### ~ USER INPUT HERE ~
edge_guess = 136770.0#117575.0#  # center of the range of ring radii
edge_lower = 136730.0#117450.0#  # lower limit of radii to consider for fit
edge_upper = 136810.0#117700.0#  # upper limit of radii to consider for fit
ring_name = 'a_ring' #'b_ring'#  # ring designation to use in output file name
path = '../output/'
### ~ END USER INPUT ~
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
outform = ["%22s","%14s","%14s","%28s","%10s","%10s","%13s","%8s","%8s","%1s"]
sl = [(6-len(df))+int(df[1:-1]) for df in outform]
outform = ",".join(f.rjust(6) for f in outform)
# output file
out = open(ring_name+'_edges.csv','w')
headers = ['OBS_ID','OET (SPM)','RET (SPM)','OET (UTC)','RDOT (km/s)',
            'Long (deg)','Edge (km)','Err (km)','RMS Resid','Flag']
out.write(",".join(h.rjust(s) for h,s in zip(headers,sl))+"\n")
plotdir = '../output/'+ring_name+'/'
if not os.path.exists(plotdir):
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
            # convert to strings for output
            oet = str(round(rf.cent_oet_spm,6))                 #
            oet = ' '*(14-len(oet))+oet             # 14 chars
            ret = str(round(rf.cent_ret_spm,6))                 #
            ret = ' '*(14-len(ret))+ret             # 14 chars
            rirv  = str(round(rf.rho_dot_kms,6))              #
            rdot = ' '*(10-len(rirv))+rirv           # 9 chars
            ilong = str(round(rf.ilong_deg,6))             #
            ilong = ' '*(10-len(ilong))+ilong       # 10 chars
            edge = str(round(rf.cent_km,6))             #
            edge = ' '*(13-len(edge))+edge
            edge_err = str(round(rf.cent_km_err,6))     #
            rms = str(round(rf.rms_resid,6))
            flag = str(rf.flag)
            # tuple of values to insert into formatted string
            row = (rf.obsid,oet,ret,rf.cent_oet_utc,rirv,ilong,edge,edge_err,rms,flag)
            # write formatted string to file
            out.write(outform % row +'\n')
            # plot
            if len(rf.rho) > 0:
                plt.plot(rf.rho,rf.pow,'-k',lw=1)
                plt.axvline(rf.cent_km,color='C1')
                plt.plot(rf.rho,rf.fit,'-r')
                plt.xlabel(r'$\rho$ (km)')
                plt.ylabel(r'$1-\exp(-\tau)$')
                plt.title(rf.obsid+' '+ring_name)
                plt.tight_layout()
                plt.savefig(plotdir+rf.obsid+'_ring_fit.png')
                plt.close()
# close file
print('')
out.close()
