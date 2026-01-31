# batch_parallel_raw2DLP_1kHz_to_50m_20260130.py
#
# Raw to GEO,CAL,DLP at variety of resolutions using 1 kHz files
# dr_km_desireds= [0.5,0.25,0.1,0.05,0.025]# Radial spacing between points (or half the labeled DLP resolution
# No diffraction reconstruction performed.
# No decimation

import sys
import time
import traceback
import rss_ringoccs as rss
from rss_ringoccs.tools.history import write_history_dict
from rss_ringoccs.tools.write_output_files import write_output_files
from multiprocessing import Process
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import fcntl

# ***** Pipeline inputs *****

### Global inputs
verbose = False
write_file_GEO = True
write_file_CAL = True
write_file_DLP = True

### RSRReader
rsr_file_list = '../tables/rsr_1kHz_files_before_USO_failure.txt' # note that there are a few 16 kHz files where no 1 kHz exists

local_path_to_output = '../output/'
local_path_to_tables = '../tables/'
program = 'batch_parallel_raw_to_DLP_1kHz_to_50m_20260130'

logfile = local_path_to_tables + program +'_log.txt' 

decimate_16khz_to_1khz = False       # Decimate 16 kHz rsr file to 1 kHz
with16 = True # load a 16kHz file if a 1 kHz file is not available

### Geometry
kernels = '../tables/e2e_kernels.ker'  # Path to meta-kernel or list of paths to
                                       #       necessary kernel files
planet = 'Saturn'                      # Name of target planet
spacecraft = 'Cassini'                 # Name of target spacecraft
pt_per_sec = 1.0                       # Number of points per second calculated
                                       #       for each parameter
### Calibration
dt_cal = 1.0                        # Time spacing in seconds between points
pnf_order = 3                       # Power normalization fit order
interact = False                    # Manually update power normalization fit

### DiffractionLimitedProfile

dr_km_desireds= [0.5,0.25,0.1,0.05,0.025]# Radial spacing between points (or
#                                    #       half of desired DLP resolution)

profile_range = [65000., 150000.]   # Lower and upper radial bounds to the
                                    #       output diffraction profile
# ***************************

contents = ['../data/'+line.strip('\n') for line in open(
         rsr_file_list,'r').readlines()]

# strip lines with # (comments)
files = []
for file in contents:
    if '#' not in file:
        files.append(file)

nfiles = len(files)

init_time = time.time()

def task(ind):
    print('\nn='+str(ind))
    rsr_file = files[ind]

    # exclude 16 kHz files? except if with16=True
    if rsr_file[-1] == '2' and not with16:
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        #continue
        return

    try:
        st = time.time()

        # print RSR file
        print(rsr_file)

        # Create instance with rsr file contents
        rsr_inst = rss.rsr_reader.RSRReader(
            rsr_file,
            verbose = verbose,
            decimate_16khz_to_1khz = decimate_16khz_to_1khz
        )
        

        # Create instance with geometry parameters
        geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft,
                kernels,local_path_to_tables=local_path_to_tables,
                local_path_to_output=local_path_to_output,
                verbose=verbose, write_file=write_file_GEO)
      
        geo_file= geo_inst.outfiles[0]+'.TAB'

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
                    verbose=verbose, write_file=write_file_CAL, 
                    local_path_to_output=local_path_to_output,
                    pnf_order=pnf_order, interact=False)
        cal_file= cal_inst.outfiles[0]+'.TAB'

        # Create instance with diffraction-limited profile and other
        #   inputs needed for diffraction correction, looping over requested resolutions
        print('logfile',logfile)
        f = open(logfile,"a")
        f.write(geo_file + '\n' + cal_file + '\n')
        f.close()
        
        for dr_km_desired in dr_km_desireds:
            dlp_inst_ing, dlp_inst_egr = (
                rss.calibration.DiffractionLimitedProfile.create_dlps(
                    rsr_inst, geo_inst, cal_inst, dr_km_desired,
                    profile_range = profile_range,
                    local_path_to_output=local_path_to_output,
                    write_file=write_file_DLP, verbose=True))
            f = open(logfile,"a")
            if type(dlp_inst_ing) != type(None):
                dlp_ing_file = dlp_inst_ing.outfiles[0]+'.TAB'
                f.write(dlp_ing_file+ '\n')
                dr_km_returned = dlp_inst_ing.dr_km
                if verbose:
                    print('DLP dr_km_desired,dlp_inst_ing.dr_km',dr_km_desired,dlp_inst_ing.dr_km)

            if type(dlp_inst_egr) != type(None):
                dlp_egr_file = dlp_inst_egr.outfiles[0]+'.TAB'
                f.write(dlp_egr_file+ '\n')
                dr_km_returned = dlp_inst_egr.dr_km
                if verbose:
                    print('DLP dr_km_desired,dlp_inst_egr.dr_km',dr_km_desired,dlp_inst_egr.dr_km)
            f.close()  
            if dr_km_returned !=  dr_km_desired:
                print('DLP dr_km_desired,dr_km_returned',dr_km_desired,dr_km_returned)
                print('skipping remaining requested dr_km_desireds')
                break # assumes dr_km_desireds in decreasing order
        et = time.time()
        run_time = (et-st)/60.
        print('File processing time (min): ' + f'{run_time:0.2f}')
    except KeyboardInterrupt:
        sys.exit()
    except:
        tb = traceback.format_exc()
        err_file = 'pipeline_test' + time.strftime("_%Y%m%d-%H%M%S") + '.err'
        fail_file = open('../output/' + err_file, 'w')
        print("%d: Failed. Printing error message to ../output/%s"
              % (ind, err_file))
        fail_file.write('-'*48+'\n')
        fail_file.write('  '+rsr_file+'\n'+'-'*36+'\n\n'+ 'n='+
                         str(ind)+'\n'+ tb+'\n')
        fail_file.close()
# end of task
def main():
    # start the process pool
    # reduce the number so that memory is not an issue, even on Loon with 128 GB

    max_workers = 5 # on Loon - still get computer crash, presumably from memory overflow

    with ProcessPoolExecutor(max_workers = max_workers) as executor:
        # submit many tasks
        futures = [executor.submit(task,arg) for arg in range(nfiles)]
        print('Waiting for tasks to complete...')
        # update each time a task finishes
        for arg in as_completed(futures):
            # report the number of remaining tasks
            print(f'About {len(executor._pending_work_items)} tasks remain')

    print('Done', flush=True)

    final_time = time.time()
    total_batch_time = (final_time - init_time)/60.
    print('Total processing time (minutes): ' + f'{total_batch_time:0.2f}')
if __name__ == '__main__':
    main()
