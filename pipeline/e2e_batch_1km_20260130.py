# This program should be run using the following command:
# echo y | python3 e2e_batch_1km_20260130.py
# This is also contained in the shell script python3 e2e_batch_1km_20260130.sh,
# which may be executed simply by typing e2e_batch_1km_20260130.sh

program = 'e2e_batch_1km_20260130'

import rss_ringoccs_global_paths # required to allow program to run from other directories

global global_path_to_rss_ringoccs
global_path_to_rss_ringoccs = rss_ringoccs_global_paths.global_path_to_rss_ringoccs

global global_path_to_local # note shortened name in next line
global_path_to_local = rss_ringoccs_global_paths.global_path_to_rss_ringoccs_local

# auxiliary path definitions based on above two globals

global global_path_to_data
global_path_to_data = global_path_to_rss_ringoccs + 'data/'

global global_path_to_output
global_path_to_output = global_path_to_rss_ringoccs + 'output/'

global global_path_to_tables
global_path_to_tables = global_path_to_rss_ringoccs + 'tables/'

global global_path_to_kernels
global_path_to_kernels = global_path_to_rss_ringoccs + 'kernels/'

global global_path_to_demo
global_path_to_demo = global_path_to_rss_ringoccs + 'demo/'

global global_path_to_demo_figs
global_path_to_demo_figs = global_path_to_rss_ringoccs + 'demo/figs/'

global global_path_to_local_output
global_path_to_local_output = global_path_to_local + 'output/'

global global_path_to_local_data
global_path_to_local_data = global_path_to_local + 'data/'

global global_path_to_local_figs
global_path_to_local_figs = global_path_to_local + 'program/figs/'

global global_path_to_local_picklefiles
global_path_to_local_picklefiles = global_path_to_local + 'picklefiles/'

global global_path_to_local_tables
global_path_to_local_tables = global_path_to_local + 'tables/'

global global_path_to_local_tmp
global_path_to_local_tmp = global_path_to_local + 'tmp/'

import sys
import time
import traceback
import rss_ringoccs as rss
from rss_ringoccs.tools.history import write_history_dict
from rss_ringoccs.tools.write_output_files import write_output_files
from rss_ringoccs.tools.pds3_tau_series import get_rev_info_from_dlp

from rss_ringoccs_local_tools import *

local_path_to_output = global_path_to_output # '../output/'
local_path_to_tables = global_path_to_tables #'../tables/'
local_path_to_data = global_path_to_data

# ***** Pipeline inputs *****
### Global inputs
verbose = False
write_file = True
### RSRReader
rsr_file_list = local_path_to_tables+'rsr_1kHz_files_before_USO_failure.txt'
decimate_16khz_to_1khz = False       # Decimate 16 kHz rsr file to 1 kHz
with16 = True				# if there is a 16 kHz file in the list, go ahead and use it

### Geometry
kernels = local_path_to_tables+'e2e_kernels.ker'  # Path to meta-kernel or list of paths to
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
dr_km_desired = 0.35 #              # Radial spacing between points (or
                                    #       half of desired DLP resolution)
profile_range = [65000., 150000.]   # Lower and upper radial bounds to the
                                    #       output diffraction profile
### DiffractionCorrection
res_km =1.0 				# Reconstruction resolution
resolution_factor = 0.75            # Factor to be multiplied to resolution
                                    #       to follow MTR86 definition


inversion_range = [74000,145000] #profile_range     # standard range for EAM PDS runs
psitype = 'fresnel' #'fresnel' #'newtonfilon11' #'fresnel' # newtonfilon11'                # Psi type
wtype = 'kbmd20'                    # Window type
sigma = 2.e-13                      # Allan deviation
fwd = False                         # Compute forward model
norm = True                         # Normalize reconstructed complex
                                    #       transmittance by window width
bfac = False                        # Use input sigma in window calculation

sat_radius =  60268.
rings_km = [74490., 91983., 117516., 122052., 136774., 139826.]
radii = [74.490, 91.983, 117.516, 122.052, 133.424, 136.774, 140.461]
Saturn = 699
Earth   = 399
Cassini = -82
ref = 'J2000'
planet = Saturn


Rvals=[74490., 91983.,117516., 122052.,136774.]
lw1 = 1.0
# ***************************

# Create new error file, removed at end if empty
err_file = program+ time.strftime("_%Y%m%d-%H%M%S") + '.err'
err_file_path = local_path_to_output + err_file
fail_file = open(err_file_path, 'w')

files_ = [local_path_to_data+line.strip('\n') for line in open(
         rsr_file_list,'r').readlines()]
files = []
for file in files_:
    if '#' in file:
        continue
    else:
        files.append(file)
 

nfiles = len(files)
init_time = time.time()

for ind in range(nfiles):
    print(ind,files[ind])
# to specify a specific rsr file:
#    if ind <= 60:
#        continue

    print('\nn='+str(ind))
    rsr_file = files[ind]

    # exclude 16 kHz files (but setting with16=True permits Rev008 16 kHz file included in 1kHz rsr file list
    if rsr_file[-1] == '2' and not with16:
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        continue

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
        geo_inst = rss.occgeo.Geometry(
            rsr_inst,
            'Saturn',
            spacecraft,
            kernels,
            local_path_to_tables=local_path_to_tables,
            local_path_to_output=local_path_to_output,
            verbose = verbose,
            write_file = write_file
        )
        geo_file = geo_inst.outfiles[0]+'.TAB'

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(
            rsr_inst,
            geo_inst,
            verbose = verbose,
            write_file = write_file,
            local_path_to_output=local_path_to_output,
            pnf_order = pnf_order,
            interact=interact
        )
        cal_file = cal_inst.outfiles[0]+'.TAB'
        # Create instance with diffraction-limited profile and other
        #   inputs needed for diffraction correction
        dlp_inst_ing, dlp_inst_egr = (
            rss.calibration.DiffractionLimitedProfile.create_dlps(
                rsr_inst,
                geo_inst,
                cal_inst,
                dr_km_desired,
                profile_range = profile_range,
                write_file = write_file,
                local_path_to_output=local_path_to_output,
                verbose = verbose
            )
        )

        # Invert profile for full occultation
        if dlp_inst_ing is not None:
            tstart = time.time()
            dlp_inst = dlp_inst_ing
            data = data_from_inst(dlp_inst)
            tau_inst = rss.DiffractionCorrection(
                data,
                res_km,
                rng = inversion_range,
                resolution_factor = resolution_factor,
                psitype = psitype,
                wtype = wtype,
                use_fwd = fwd,
                use_norm = norm,
                bfac = bfac,
                verbose = verbose
            )
            tend = time.time()

            tau_history = set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,resolution_factor,psitype,wtype,program,rssocc_version='1.3-beta')
            tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
            if write_file:
                rev_info = dlp_inst.rev_info
                tau_files = write_output_files.write_output_files(tau_inst,rev_info=rev_info,
                history = tau_history, local_path_to_output = local_path_to_output)

                taufile_PDS = get_CORSS_8001_TAUfile(rev_info)
                print('Creating summary PDF...')
                rss.tools.plot_summary_doc_v5(geo_inst, cal_inst, dlp_inst, tau_inst, tau_files,
                    psitype,res_km,taufile_PDS=taufile_PDS,
                    reslocs_sav = global_path_to_local_tables + 'resloc_v6.sav',
                    wavefile = global_path_to_local_tables + 'wave_list_20260129.csv'
                )

        if dlp_inst_egr is not None:
            tstart = time.time()
            dlp_inst = dlp_inst_egr
            data = data_from_inst(dlp_inst)
            tau_inst = rss.DiffractionCorrection(
                data,
                res_km,
                rng = inversion_range,
                resolution_factor = resolution_factor,
                psitype = psitype,
                wtype = wtype,
                use_norm = norm,
                bfac = bfac,
                verbose = verbose
            )
            tend = time.time()

            # tau_inst.history = write_history_dict(tau_inst.input_vars,
            #                                       tau_inst.input_kwds,
            #                                       __file__)
            tau_history = set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,resolution_factor,psitype,wtype,program,rssocc_version='1.3-beta')
            tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
            if write_file:
                rev_info = dlp_inst.rev_info
                tau_files = write_output_files.write_output_files(tau_inst,rev_info=rev_info,
                history = tau_history, local_path_to_output = local_path_to_output)

                taufile_PDS = get_CORSS_8001_TAUfile(rev_info)
#                rss.tools.plot_summary_doc_v4(
                print('Creating summary PDF...')
                rss.tools.plot_summary_doc_v5(geo_inst, cal_inst, dlp_inst, tau_inst, tau_files,
                    psitype,res_km,taufile_PDS=taufile_PDS,
                    reslocs_sav = global_path_to_local_tables + 'resloc_v6.sav',
                    wavefile = global_path_to_local_tables + 'wave_list_20260129.csv'
                )

        et = time.time()
        run_time = (et-st)/60.
        print('File processing time (min): ' + f'{run_time:0.2f}')
    except KeyboardInterrupt:
        sys.exit()
    except:
        tb = traceback.format_exc()
        print("%d: Failed. Printing error message to ../output/%s"
              % (ind, err_file))
        fail_file.write('-'*48+'\n')
        fail_file.write('  '+rsr_file+'\n'+'-'*36+'\n\n'+ 'n='+
                         str(ind)+'\n'+ tb+'\n')

final_time = time.time()
total_batch_time = (final_time - init_time)/60./60.
print('Total processing time (hrs): ' + f'{total_batch_time:0.3f}')
fail_file.close()
if os.path.getsize(err_file_path) == 0:
        os.remove(err_file_path)
print("All done!")
