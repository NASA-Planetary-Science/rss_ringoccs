# The $-delimited keywords in this template are filled by the driver program
# construct_batch_parallel_CSVtoTAU_20260130.ipynb
# to create a sstandalone Python program that uses multiprocessing to produce
# TAU files for the specified Rev, band, DSN, direction, radial range, psitype,
# resolution, resolution_factor over the specified inversion_range

program = 'batch_parallel_CSVtoTAU_template-20260130a'

description = '$DESCRIPTION$' # to make it easier to trace what this run was 

from rss_ringoccs_local_tools import * 

def processing_range_required(inversion_range,rkm,F,res_factor,res_km):
    if(inversion_range[0] <np.min(rkm)) or (inversion_range[1] > np.max(rkm)):
        raise Exception('GEO file does not span requested inversion range')
    W = processing_window(F,res_factor,res_km) # processing window across GEO file
    fWofrkm= interpolate.interp1d(rkm,W)
    min_proc_range_req = int(np.floor(inversion_range[0]-fWofrkm(inversion_range[0])/2))
    max_proc_range_req = int(np.ceil(inversion_range[1]+fWofrkm(inversion_range[1])/2))
    proc_range_req = [min_proc_range_req,max_proc_range_req]
    return proc_range_req 
    
def is_processing_range_valid(inversion_range,rkm,F,res_factor,res_km,verbose=False):
    proc_range_req = processing_range_required(inversion_range,rkm,F,res_factor,res_km)
    proc_range_avail = [np.min(rkm),np.max(rkm)]
    if verbose:
        print('inversion range:',inversion_range)
        print('proc_range_req:',proc_range_req)
        print('proc_range_avail:',proc_range_avail)
    return (proc_range_avail[0]<=proc_range_req[0]) and (proc_range_avail[1]>=proc_range_req[1])

def is_processing_range_valid_from_GEO_DLP_files(inversion_range,GEOfilepath,DLPfilepath,res_factor,res_km):
    rkmDLP = np.loadtxt(dlp_file,delimiter=',',usecols=[0],unpack=True)
    rkm, rdot, F = np.loadtxt(GEOfilepath,delimiter=',',usecols=[3,8,10],unpack=True)
    if 'I_GEO' in geo_file:
        direc = 'I'
        L = np.where((rdot<0) & (rkm >= np.min(rkmDLP)) & (rkm <= np.max(rkmDLP)))[0]
    elif 'E_GEO' in geo_file: 
        direc = 'E'
        L = np.where((rdot>0) & (rkm >= np.min(rkmDLP)) & (rkm <= np.max(rkmDLP)))[0]    
    else:
        raise Exception("Unable to parse direction from GEOfilepath")
    return is_processing_range_valid(inversion_range,rkm[L],F[L],res_factor,res_km)

def processing_window(F,res_factor,res_km):
    W = 2*F**2/(res_factor * res_km)
    return W

def update_inversion_range(inversion_range,geo_file,dlp_file,res_factor,res_km):
    rkmDLP = np.loadtxt(dlp_file,delimiter=',',usecols=[0],unpack=True)
    rkm, rdot, F = np.loadtxt(geo_file,delimiter=',',usecols=[3,8,10],unpack=True)
    if 'I_GEO' in geo_file:
        direc = 'I'
        L = np.where((rdot<0) & (rkm >= np.min(rkmDLP)) & (rkm <= np.max(rkmDLP)))[0]
    elif 'E_GEO' in geo_file: 
        direc = 'E'
        L = np.where((rdot>0) & (rkm >= np.min(rkmDLP)) & (rkm <= np.max(rkmDLP)))[0]
    else:
        raise Exception("Unable to parse direction from GEOfilepath")
    inversion_range_updated = np.zeros(2,dtype=int)
    rkmL = rkm[L]
    
    FL = F[L]
#    print('FL',FL)
    # plt.plot(FL)
    WL = processing_window(FL,res_factor,res_km)
    rkmL_minus_WL = rkmL-WL/2 # at each radial location, compute corresponding min radius required for this window
    rkmL_plus_WL  = rkmL+WL/2 # at each radial location, compute corresponding min radius required for this window
    min_rkmL = np.min(rkmL) # minimum radius available
    max_rkmL = np.max(rkmL) # maximum radius available
    LL = np.where((rkmL_minus_WL >= np.min(rkmL)) & (rkmL_plus_WL <= np.max(rkmL)))[0]
    inversion_range_updated[0] = 0+np.max([inversion_range[0],np.ceil(np. min(rkmL[LL])).astype(int)])
    inversion_range_updated[1] = 0+np.min([inversion_range[1],np.floor(np.max(rkmL[LL])).astype(int)])
    if inversion_range_updated[1]<=inversion_range_updated[0]:
        print('Illegal update to',inversion_range,':',inversion_range_updated)
        raise Exception('No part of requested inversion range can be processed for this event')
    return inversion_range_updated
    
init_time = time.time()

$RES_FACTOR$ # res_factor = 0.75 # so resolution agrees with PDS
wtype = '$WTYPE$' # 'kbmd20' # modify this only for testing of window functions

tmp2output = True # copy files to rss_ringoccs_local/output/
tmp2globaloutput = True # instead, send to rss_ringoccs/output/
add_inversion_range = True
clean = True # move files to old/
$SILENT$
psitypes = ['$PSITYPE$']
'$RES_FACTOR$'
$RES_KM$ # res_km = value
min_dlp_res = min_required_dlp_res_km(res_km,res_factor)*1000 # meters

# There is a standard set of DLP files at the following resolutions
# Choose the one with the coarsest allowable resolution to speed up
# execution 

if min_dlp_res <20:
    dlp_res='DLP_*10M_'
elif min_dlp_res <40:
    dlp_res='DLP_*20M_'
elif min_dlp_res <50:
    dlp_res='DLP_*40M_'
elif min_dlp_res <100:
    dlp_res='DLP_*50M_' 
elif min_dlp_res <200:
    dlp_res='DLP_*100M_'
elif min_dlp_res <500:
    dlp_res='DLP_*200M_'
else:
    dlp_res='DLP_*500M_'

output_creation_date = '$OUTPUT_CREATION_DATE$' # ex: '2025122*' # Date of DLP file

search_string = output_creation_date+'_0001.TAB' 
alt_search_string = output_creation_date+'_0002.TAB'

rev = '$REV$'
direc = '$DIREC$'
band = '$BAND$'
dsn = '$DSN$'
name = '$NAME$'
$DRHO_KM$

$RES_FACTOR$ #line= line.replace('$RES_FACTOR$','res_factor = '+str(res_factor))

# the processing_range is the range of the DLP file to use, NOT the inversion_range
processing_range = [int('$PROCESSING_RANGE_MIN$'),int('$PROCESSING_RANGE_MAX$')]

SEP = os.sep

search_dir = global_path_to_output + SEP + 'Rev'+rev+ SEP+'Rev'+rev+'*'+direc+SEP+'*'
dirs = glob.glob(search_dir)

results_all = []
#bands = np.array(run['bands'])

string = band + dsn+'_'+direc
# find the correct output directory that contains GEO/TAB/DLP for this rev, direction, band, dsn 
for this_dir in dirs:
    if string in this_dir:
        # if not silent:
        #     print('Processing',this_dir)
        if not silent:
            print('now search in ',this_dir,'for GEO,CAL,DLP')
        GEO_glob = SEP+'*GEO*'+search_string
        CAL_glob =  SEP+'*CAL*'+search_string
        DLP_glob =  SEP+'*'+dlp_res+search_string        
        try:
            geo_file = glob.glob(this_dir + GEO_glob)[0]
        except:
            raise Exception('No matching GEO file found for '+GEO_glob)
        try:
            cal_file = glob.glob(this_dir + CAL_glob)[0]
        except:
            raise Exception('No matching CAL file found for '+CAL_glob)
        try:
            dlp_file = glob.glob(this_dir + DLP_glob)[0]
        except:
            print('No matching DLP file at requested resolution',dlp_res)
            print(this_dir + DLP_glob)
            dlp_res = 'DLP_*???M'
            DLP_glob_all = SEP+'*'+dlp_res+'_'+alt_search_string
            dlp_files = glob.glob(this_dir + DLP_glob_all)
            if verbose:
                print('found:')
                for dlp_file in dlp_files:
                    print(dlp_file)
            try:
                dlp_file = dlp_files[0] # use highest resolution one in list
            except:
                raise Exception('No matching DLP file found')
        #print('Using',dlp_file)
        if not silent:
            print(geo_file)
            print(cal_file)
            print(dlp_file)

        # update the inversion range based on the available range in the DLP file and the
        # required window sizes for the resolution res_km of this run
        
        inversion_range_orig = [int('$INVERSION_RANGE_MIN$'),int('$INVERSION_RANGE_MAX$')]
        inversion_range = inversion_range_orig

        try:
            inversion_range =update_inversion_range(inversion_range_orig,
                                                    geo_file,dlp_file,res_factor,res_km)
        except:
            raise Exception('Error in update_inversion_range()')
        try:
            is_processing_range_valid=is_processing_range_valid_from_GEO_DLP_files(inversion_range,
                                                            geo_file,dlp_file,res_factor,res_km)
        except:
            raise Exception('Error in is_processing_range_valid_from_GEOfile()')

        if is_processing_range_valid==False:
            raise Exception("processing range is invalid")
        if not silent and (inversion_range_orig !=  inversion_range.tolist()):
            print('*** Inversion range modified from ',\
                  inversion_range_orig,'to',inversion_range.tolist())

# get radial range of dlp_file to see if reqested range is present
        dlp_contents = np.loadtxt(dlp_file,delimiter=',')
        rmin_dlp = dlp_contents[0,0]
        rmax_dlp = dlp_contents[-1,0]
        if (rmin_dlp > inversion_range[0]) or (rmax_dlp < inversion_range[0]):
            print('DLP range does not span requested inversion_range')
            print('DLP range:',rmin_dlp,rmax_dlp)
            print('inversion_range:',inversion_range)
            raise Exception('Error due to inadequate DLP range')
        if not silent: # this gets printed for each CPU!
            print(os.path.basename(geo_file))
            print(os.path.basename(cal_file))
            print(os.path.basename(dlp_file))
            
        $TRIM_DLP$ # trim_dlp = True or False
        if trim_dlp:
            dlp_file= trim_dlp_file(dlp_file,processing_range,
                                    verbose=not silent,write=True,overwrite=False)
        title = name + ' Rev' + rev + direc 
        this_start_time = time.time()
        data = rss_ringoccs.ExtractCSVData(geo_file, cal_file, dlp_file)
        this_stop_time = time.time()
        if drho_km != 0:
            data.rho_km_vals += drho_km
        psitype = psitypes[0]
        
$MAX_WORKERS$ # max_workers = cpu_count()

#drange is the chunck size in km of the inversion_ranged processed by each cpu

drange = math.ceil((inversion_range[1]-inversion_range[0])/max_workers)

nranges,min_ranges,max_ranges = subdivide_inversion_range(inversion_range, drange)

def task(ind):
    include_history=False
    try:
        this_inversion_range = [int(min_ranges[ind]),int(max_ranges[ind])]
        this_start_time = time.time()
        
        tau_inst = rss_ringoccs.DiffractionCorrection(
               data, res_km, rng=this_inversion_range, resolution_factor=res_factor,
               psitype=psitype, wtype=wtype, verbose=False)
        
        tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
        
        this_stop_time = time.time()
        
        print("Diffraction-correction processing time for this chunk:",
             this_stop_time - this_start_time,"seconds")
        
        rev_info = get_rev_info_from_dlp(dlp_file)

        if include_history:
             tau_history = set_tau_history_from_csv(tau_inst,geo_file,cal_file,dlp_file,
                  this_start_time,this_stop_time,res_km,
                  this_inversion_range,res_factor,psitype,wtype,program,
                  rssocc_version='1.3-beta')
        else: # specify a legal empty dictionary
                    tau_history = {
            'key_order0': ['User Name', 'Host Name', 'Operating System',
                        'Python Version', 'rss_ringoccs Version']
            ,'key_order1': ['Source Directory','Source File', 
                        'Positional Args', 'Keyword Args', 'Additional Info']
            , 'hist name': 'DiffractionReconstruction history'
            , 'User Name': ''
            , 'Host Name': ''
            , 'Run Date': ''
            , 'Python Version': ''
            , 'rss_ringoccs Version': ''
            , 'Operating System': ''
            , 'Source Directory': ''
            , 'Source File': ''
            , 'Positional Args': ''
            , 'Keyword Args': ''
            , 'Additional Info': ''
            , 'description': ''
            }           

        outfiles = write_output_files.write_output_files(tau_inst,rev_info=rev_info,
                add_suffix = f'.{ind+1:04d}',# add index+1 to filename so that it is unique
                history = tau_history, local_path_to_output = global_path_to_local_tmp)
    except:
        tb = traceback.format_exc()
        err_file = program + time.strftime("_%Y%m%d-%H%M%S") + '.err'
        err_filepath = global_path_to_output + err_file

        fail_file = open(err_filepath, 'w')
        print("%d: Failed. Printing error message to %s" % (ind, err_filepath))
        fail_file.write('-'*48+'\n')
        fail_file.write('  '+dlp_file+'\n'+'-'*36+'\n\n'+ 'n='+
                     str(ind)+'\n'+ tb+'\n')
        fail_file.close()
# end of task

def main():
# pre-clean the temporary directory (if it exists!) of any previous TAU files so they won't get merged with current run
    if clean:
        search_dir = global_path_to_local_tmp +'/Rev'+rev+'/Rev'+rev+'*'+direc+'/'+'*'+band+dsn+'*'+direc+'/'
        indirs = glob.glob(search_dir)
        if len(indirs) != 0:
            indir = indirs[0]
            try:
                result = subprocess.run("ls "+ indir+"RSS*TAU*", shell=True,capture_output=True, text=True, check=True)
                output_lines = result.stdout.strip().splitlines()
                for line in output_lines:
                    print(os.path.basename(line))
                    cmd = "mv "+line+" "+indir+"old/"
                    print(cmd)
                    result = subprocess.run(cmd, shell=True,capture_output=True, text=True, check=True)
            except:
                print('No stale TAU files found in '+indir)
        else:
            print("Temporary directory for output files does not (yet) exist")

#    start the process pool
    os.system('date')
    $MAX_WORKERS$ # max_workers = cpu_count()
    with ProcessPoolExecutor(max_workers = max_workers) as executor:
#       submit many tasks
         futures = [executor.submit(task,arg) for arg in range(nranges)]
         print('Waiting for '+str(nranges)+' multiprocessor tasks to complete...')
#        update each time a task finishes
         for arg in as_completed(futures):
#            report the number of remaining tasks
             print(f'About {len(executor._pending_work_items)} tasks remain')
#       submit many tasks  
    search_dir = global_path_to_local_tmp +'/Rev'+rev+'/Rev'+rev+'*'+direc+'/'+'*'+band+dsn+'*'+direc+'/'
    #print('\nin main(): search_dir = ',search_dir)
    indir = glob.glob(search_dir)[0]
    outfile,outlblfile,outfile_psitype,outlblfile_psitype,outdir = \
        merge_sorted_tabfiles(indir,inversion_range,tmp2output=tmp2output,psitype=psitype,
                              tmp2globaloutput=tmp2globaloutput,update_radius_corrections = False,
                              add_inversion_range=add_inversion_range,clean=clean,silent=False)
    print('\nMerged TAU files saved in \n',outdir,flush=True)
    for file in (outfile,outlblfile,outfile_psitype,outlblfile_psitype):
        print(os.path.basename(file),flush=True)
    
    print('All Done!', flush=True)
    os.system('date')
    
    final_time = time.time()
    total_batch_time = (final_time - init_time)/60.
    print('Total processing time: ' + f'{total_batch_time:0.2f} minutes')
    print(get_processor_info())

if __name__ == '__main__':
    main()