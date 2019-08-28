### ~ USER INPUT HERE ~
edge_guess = 136770.0#117575.0#  # center of the range of ring radii
edge_lower = 136730.0#117450.0#  # lower limit of radii to consider for fit
edge_upper = 136810.0#117700.0#  # upper limit of radii to consider for fit
ring_name = 'a_ring_postUSO' #'b_ring'#  # ring designation to use in output file name
path = '../output/'#postUSO_jolene/'
#path = '/Volumes/sflury001/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/'
plot = True
### ~ END USER INPUT ~
import numpy as np
import sys, os
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import matplotlib.pyplot as plt
# use the os.listdir package to find all directories in rss_ringoccs
# output file directory structure
occdirs = [path+dir1+'/'+dir2+'/'+dir3+'/'
            for dir1 in os.listdir(path)
            if '.' not in dir1
            for dir2 in os.listdir(path+dir1)
            if '.' not in dir1+'/'+dir2
            for dir3 in os.listdir(path+dir1+'/'+dir2)
            if '.' not in dir1+'/'+dir2+'/'+dir3]
#
# output format:
#   Observation ID, observed event time in SPM of ring edge, UTC at
#   which edge was observed, resolution of profile, intertial longitude
#   at which ring edge was observed, ring event time at which edge
#   received spacecraft signal, best fit ring edge (km), uncertainty in
#   ring edge, fit result, and the sum squared residuals
#
outform = ["%22s","%14s","%14s","%28s","%12s","%12s","%13s","%8s","%8s","%12s","%1s"]
sl = [(6-len(df))+int(df[1:-1]) for df in outform]
outform = ",".join(f.rjust(6) for f in outform)
# output file
out = open('../tables/'+ring_name+'_edges.csv','w')
#rfo = open('KLTS.data.combined.RSS.cassini.052','w')
headers = ['OBS_ID','OET (SPM)','RET (SPM)','OET (UTC)','RDOT (km/s)',
            'Long (deg)','Edge (km)','Err (km)','RMS Resid','R drift (km)','Flag']
out.write(",".join(h.rjust(s) for h,s in zip(headers,sl))+"\n")
# directory for plots
plotdir = '../output/'+ring_name+'/'
if not os.path.exists(plotdir) and plot:
    os.system('[ ! -d ' + plotdir + ' ] && mkdir -p ' + plotdir)
# iterate over all directories
for n,dir in enumerate(occdirs):
    # iterate over all TAU files in the directory
    for file in os.listdir(dir):
        if file.endswith('.TAB') and 'TAU' in file:
            # occ name
            occ = dir.split('/')[-2].split('_')[0]
            # skip occultations not published on PDS
            if occ == 'Rev053CE' :
                continue

            # call rss_ringoccs ring fitter
            rf = rss.tools.ring_fit(dir+file,edge_guess,data_lims=[edge_lower,edge_upper])
            # find and input geo file
            splts = np.array([splt for splt in file.split('_')],dtype=str)
            geo_file = "_".join(splt for splt in splts[:5])+'_GEO_'+"_".join(splt for splt in splts[7:])
            if not os.path.exists(dir+geo_file):
                if splts[-1][-5] == '1':
                    splts[-2] = str(int(splts[-2])-1)
                else:
                    splts[-1] = '000'+str(int(splts[-1][-5])-1)+'.TAB'
                geo_file = "_".join(splt for splt in splts[:5])+'_GEO_'+"_".join(splt for splt in splts[7:])

            goet,rho_dot,Fresnel = np.loadtxt(dir+geo_file,delimiter=',',usecols=(0,8,10)).T
            # extract ring intercept velocity by interpolating at ring edge
            rho_dot_kms = np.interp(rf.cent_oet_spm,goet,rho_dot)
            # obtain radial offset due to phase drift
            if len(rf.rho[(rf.rho>rf.cent_km)]) > 3:
                # obtain and interpolate Fresnel scale to edge fit
                F = np.interp(rf.cent_oet_spm,goet,Fresnel)
                phase = np.unwrap(rf.phase*np.pi/180.)
                p = np.polyfit(rf.rho[(rf.rho>rf.cent_km)],phase[(rf.rho>rf.cent_km)],1)
                # get phase-induced radial offset
                r_off = p[0]*(F**2.)/np.pi
            else:
                r_off = -9999.999999

            # check edge width
            if rf.flag == 0 :
                if abs(rf.fit_parameters[2]) < 2:
                    rf.flag = int(4)
            # convert to strings for output
            oet = str(round(rf.cent_oet_spm,6))                 #
            oet = ' '*(14-len(oet))+oet             # 14 chars
            ret = str(round(rf.cent_ret_spm,6))                 #
            ret = ' '*(14-len(ret))+ret             # 14 chars
            rirv  = str(round(rho_dot_kms,6))              #
            rdot = ' '*(12-len(rirv))+rirv           # 9 chars
            ilong = str(round(rf.ilong_deg,6))             #
            ilong = ' '*(12-len(ilong))+ilong       # 10 chars
            edge = str(round(rf.cent_km,6))             #
            edge = ' '*(13-len(edge))+edge
            edge_err = str(round(rf.cent_km_err,6))     #
            rms = str(round(rf.rms_resid,6))
            flag = str(rf.flag)
            r_off = str(round(r_off,6))
            r_off = ' '*(12-len(r_off))+r_off
            # tuple of values to insert into formatted string
            row = (rf.obsid,oet,ret,rf.cent_oet_utc,rirv,ilong,edge,edge_err,rms,r_off,flag)
            # write formatted string to file
            out.write(outform % row +'\n')

            # plot
            if len(rf.rho) > 0 and plot:
                plt.plot(rf.rho,rf.pow,'-k',lw=1)
                plt.axvline(rf.cent_km,color='C1')
                plt.axvline(rf.cent_km-float(r_off.strip('')),color='C2')
                plt.plot(rf.rho,rf.fit,'-r')
                plt.xlim(round(rf.cent_km)-50,round(rf.cent_km)+50)
                plt.xlabel(r'$\rho$ (km)')
                plt.ylabel(r'$1-\exp(-\tau)$')
                plt.title(rf.obsid+' '+ring_name)
                plt.tight_layout()
                plt.savefig(plotdir+rf.obsid+'_ring_fit_phase.png')
                plt.close()
# close file
print('')
out.close()
