#!/usr/bin/env python
# coding: utf-8

# # construct_batch_parallel_CSVtoTAU_20260130.ipynb
# 
# Sample code for pipeline to demonstrate parallel processing
#     
# 

# In[1]:


# all packages and utility programs loaded here

from rss_ringoccs_local_tools import *

# this also defines global variables 


# In[16]:


execute_this_cell = True
if execute_this_cell:

    from subprocess import DEVNULL, STDOUT

    description = 'Example code to produce TAU files for C + inner B ring for all bands, stations'

    version = '20260130a' # prepend this name to output *.sh file that contains sequence of python commands to run
    template_file = 'batch_parallel_CSVtoTAU_template-20260130a.py'
    rundir = 'runs_'+version # create subdirectory to hold template-filled *.py files for this run

    if not os.path.isdir(rundir):
        # print('The pipeline/'+rundir+'/ directory needs to exist and to contain the following symbolic links:')
        # print('rss_ringoccs_global_paths.py -> ../../rss_ringoccs_global_paths.py')
        # print('rss_ringoccs_local_tools.py -> ../../rss_ringoccs/tools/rss_ringoccs_local_tools.py')
        # print('Attempting to set up the '+rundir+'/ directory:')
        try:
            os.makedirs(rundir, exist_ok=True)
            directory_path = os.getcwd()+'/'+newdir+'/'
            command = ["ln", "-s","../../rss_ringoccs_global_paths.py"] 
            result = subprocess.run(command, cwd=directory_path, stdout=DEVNULL, stderr=DEVNULL)
            command = ["ln", "-s",".../../rss_ringoccs/tools/rss_ringoccs_local_tools.py"] 
            result = subprocess.run(command, cwd=directory_path, stdout=DEVNULL, stderr=DEVNULL)
            if not os.path.islink(directory_path +'rss_ringoccs_global_paths.py') or \
             not os.path.islink(directory_path +'rss_ringoccs_local_tools.py'):
                raise Exception("Failed to set up symbolic links")
            print('Successfully set up',rundir+'/')
        except:
            raise Exception("Unable to set up",rundir,"and required symbolic links.")

    outfile_root = rundir+'/batch_parallel_CSVtoTAU-'+version+'-'
    shell_file = 'batch_parallel_CSVtoTAU-'+version+'.sh'
    f = open(template_file,'r')
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()

    # optionally, code can produce TAU file for a given trio of GEO, CAL, DLP
    # Ex:
    # data_dir = global_path_to_output + 'Rev133/Rev133E/Rev133E_RSS_2010_170_X43_E/' # note / at end
    # geo_file = data_dir + "RSS_2010_170_X43_E_GEO_20250408_0001.TAB"
    # cal_file = data_dir + "RSS_2010_170_X43_E_CAL_20250408_0001.TAB"
    # dlp_file = data_dir + "RSS_2010_170_X43_E_DLP_0050M_20250408_0001.TAB"


    # normal case is NOT to specify GEO,CAL,DLP, so use these definitions:

    geo_file = 'NONE'
    cal_file = 'NONE'
    dlp_file = 'NONE'

    res_kms = np.array([1.0, 0.500, 0.400, 0.300, 0.250, 0.200, 0.150, 0.100, 0.075, 0.050])
    min_res_km_Sband = 0.200 # don't try to process lower resolution S band than this value

    bands = ['X','K','S']

# not all of these include the C ring but code checks for this
    revs = [
            'Rev007',
            'Rev008',
            'Rev009',
            'Rev010',
            'Rev011',
            'Rev012',
            'Rev013',
            'Rev014',
            #'Rev028', #BAD EVENT - very distant
            'Rev044',
            'Rev046',
            'Rev053',
            'Rev054',
            'Rev056',
            'Rev057',
            'Rev058',
            'Rev060',
            'Rev063',
            'Rev064',
            'Rev067',
            'Rev079',
            'Rev084',
            'Rev089',
            'Rev123',
            'Rev125',
            'Rev133',
            'Rev137'
           ]            
    # specify restrictions on DSNs to include
    include_all_DSNs = True
    # otherwise, do something like this to get 70 m stations for X band
    DSN_include=['14','43','63'] #70m stations only

    dirs = ['I','E']
    psitypes = ['NewtonFilon12']

    names = ['C + inner B ring']
    inversion_range_mins = ['74400']
    inversion_range_maxs = ['97000'] # extend a bit to get wave in inner B ring

    trim_dlp = True # to save memory in multiprocessing, trip the DLP to this range
    processing_range_mins = ['60000']
    processing_range_maxs = ['110000']
    res_factor = 0.75

    wtype = 'kbmd20'
    drho_km = 0.0 # offset applied to radius scale
    output_creation_date = '*' # if searching for specific DLP file, set date of DLP file creation

    silent = False

    index = 0 # index of first output file
    with open(shell_file,'w') as fout:
        for res_km in res_kms:
            for Rev in revs:
                for direction in dirs:
                    for band in bands:
                        if band == 'S' and res_km < min_res_km_Sband:
                            continue # S band can't handle high resolution
                        files = get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,Rev=Rev,\
                                    DSN=None,direction=direction,local_path_to_data=global_path_to_data,
                                        local_path_to_tables = global_path_to_tables,silent=True)
                        DSNs = []
                        for file in files:
                            basename = os.path.basename(file)

                            basenameband = basename[len('RSS_2005_123_K')-1]
                            basenamedir = basename[len('RSS_2005_123_K34_E')-1]
                            basenameDSN = basename[len('RSS_2005_123_K3')-1:len('RSS_2005_123_K34_')-1]

                            if not include_all_DSNs and basenameDSN not in DSN_include and band != 'K':
                                continue
                            if Rev in file and band == basenameband and direction == basenamedir:
                                DSNs.append(band+basenameDSN)
                        for DSN in DSNs:
                            for psitype in psitypes:
                                for name,inversion_range_min,inversion_range_max,processing_range_min,processing_range_max in \
                                    zip(names,inversion_range_mins,inversion_range_maxs,processing_range_mins,processing_range_maxs):
                                    f = open(template_file,'r')
                                    contents = f.read().splitlines() # strips \n at end of each line
                                    f.close()
                                    for i,line in enumerate(contents):
                                        if '$DESCRIPTION' in line:     
                                            line= line.replace('$DESCRIPTION$',description)
                                            contents[i] = line                                            
                                        if '$GEO_FILEPATH$' in line:
                                            line= line.replace('$GEO_FILEPATH$',geo_file)
                                            contents[i] = line
                                        if '$CAL_FILEPATH$' in line:
                                            line= line.replace('$CAL_FILEPATH$',cal_file)
                                            contents[i] = line
                                        if '$DLP_FILEPATH$' in line:
                                            line= line.replace('$DLP_FILEPATH$',dlp_file)
                                            contents[i] = line
                                        if '$OUTPUT_CREATION_DATE$' in line:
                                            line= line.replace('$OUTPUT_CREATION_DATE$',output_creation_date)
                                            contents[i] = line
                                        if '$WTYPE$' in line:
                                            line= line.replace('$WTYPE$',wtype)
                                            contents[i] = line
                                        if '$RES_KM$' in line:
                                            line= line.replace('$RES_KM$','res_km = '+str(res_km))
                                            contents[i] = line
                                        if '$RES_FACTOR$' in line:
                                            line= line.replace('$RES_FACTOR$','res_factor = '+str(res_factor))
                                            contents[i] = line
                                        if '$DRHO_KM$' in line:
                                            line= line.replace('$DRHO_KM$','drho_km = '+str(drho_km))
                                            contents[i] = line
                                        if '$REV$' in line:
                                            line= line.replace('$REV$',Rev[3:])
                                            contents[i] = line
                                        if '$DIREC$' in line:
                                            line= line.replace('$DIREC$',direction)
                                            contents[i] = line
                                        if '$BAND$' in line:
                                            line= line.replace('$BAND$',DSN[0])
                                            contents[i] = line
                                        if '$DSN$' in line:
                                            line= line.replace('$DSN$',DSN[1:])
                                            contents[i] = line
                                        if '$PSITYPE$' in line:
                                            line= line.replace('$PSITYPE$',psitype)
                                            contents[i] = line
                                        if '$INVERSION_RANGE_MIN$' in line:
                                            line = line.replace('$INVERSION_RANGE_MIN$',inversion_range_min)
                                            contents[i] = line
                                        if '$INVERSION_RANGE_MAX$' in line:
                                            line = line.replace('$INVERSION_RANGE_MAX$',inversion_range_max)
                                            contents[i] = line
                                        if '$NAME$' in line:
                                            line = line.replace('$NAME$',name)
                                            contents[i] = line
                                        if '$SILENT$' in line:
                                            if silent == True:
                                                line = line.replace('$SILENT$','silent = True')
                                            else:
                                                line = line.replace('$SILENT$','silent = False')
                                            contents[i] = line
                                        if '$MAX_WORKERS$' in line:
                                            line = line.replace('$MAX_WORKERS$','max_workers = cpu_count()')
                                            contents[i] = line
                                        if '$TRIM_DLP$' in line:
                                            if trim_dlp == True:
                                               line = line.replace('$TRIM_DLP$','trim_dlp = True')
                                            else:
                                               line = line.replace('$TRIM_DLP$','trim_dlp = False')
                                            contents[i] = line
                                        if '$PROCESSING_RANGE_MIN$' in line:
                                            line = line.replace('$PROCESSING_RANGE_MIN$',processing_range_min)
                                            contents[i] = line
                                        if '$PROCESSING_RANGE_MAX$' in line:
                                            line = line.replace('$PROCESSING_RANGE_MAX$',processing_range_max)
                                            contents[i] = line

                                    outfile = outfile_root + f'{index:04d}.py'
                                    print(outfile, res_km,Rev,direction,DSN,psitype,name)
                                    index += 1
                                    with open(outfile,'w') as f:
                                        for line in contents:
                                            print(line,file=f)
                                    print('python3 '+outfile,file=fout)
    print('\nExecute this shell script from the rss_ringoccs/pipeline directory:')
    print(shell_file)
    print('\nOr, to run an individual job from the rss_ringoccs/pipeline directory, type:')
    print('python3 '+outfile,'\n')
    os.chmod(shell_file,0o777)


# In[ ]:




