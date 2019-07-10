### ~ USER INPUT HERE ~
edge_guess = 136770.0#117575.#  # center of the range of ring radii
edge_lower = 136730.0#117450.#  # lower limit of radii to consider for fit
edge_upper = 136810.0#117700.#  # upper limit of radii to consider for fit
ring_name = 'a_ring_pds_test'   # ring designation to use in output file name
#path = '../output/'
path = '/Volumes/sflury001/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/'
### ~ END USER INPUT ~
import sys, os
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
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
outform = ["%22s","%14s","%14s","%27s","%10s","%10s","%13s","%8s"]
sl = [(6-len(df))+int(df[1:-1]) for df in outform]
outform = ",".join(f.rjust(6) for f in outform)
# output file
out = open(ring_name+'_edges.csv','w')
headers = ['OBS_ID','OET (SPM)','RET (SPM)','OET (UTC)','RDOT (km/s)',
            'Long (deg)','Edge (km)','Err (km)']
out.write(",".join(h.rjust(s) for h,s in zip(headers,sl))+"\n")

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
            rf = rss.tools.ring_fit(dir+file,edge_guess,edge_lims=[edge_lower,edge_upper])
            # convert to strings for output
            oet = str(round(rf.edge_oet_spm,6))                 #
            oet = ' '*(14-len(oet))+oet             # 14 chars
            ret = str(round(rf.edge_ret_spm,6))                 #
            ret = ' '*(14-len(ret))+ret             # 14 chars
            rirv  = str(round(rf.rho_dot_kms,6))              #
            rdot = ' '*(10-len(rirv))+rirv           # 9 chars
            ilong = str(round(rf.ilong_deg,6))             #
            ilong = ' '*(10-len(ilong))+ilong       # 10 chars
            edge = str(round(rf.edge_km,6))             #
            edge = ' '*(13-len(edge))+edge
            edge_err = str(round(rf.edge_km_err,6))         #
            # tuple of values to insert into formatted string
            row = (rf.obsid,oet,ret,rf.edge_oet_utc,rirv,ilong,edge,edge_err)
            # write formatted string to file
            out.write(outform % row +'\n')
            # print to terminal
            #print(rf.obsid,edge,edge_err)
# close file
print('')
out.close()
