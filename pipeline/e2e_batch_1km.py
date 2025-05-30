import sys
import time
import traceback
import rss_ringoccs as rss
from rss_ringoccs.tools.history import write_history_dict
from rss_ringoccs.tools.write_output_files import write_output_files

# ***** Pipeline inputs *****
### Global inputs
verbose = False
write_file = True
### RSRReader
rsr_file_list = '../tables/rsr_1kHz_files_before_USO_failure.txt'
decimate_16khz_to_1khz = False       # Decimate 16 kHz rsr file to 1 kHz
with16 = False

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
dr_km_desired = 0.25                # Radial spacing between points (or
                                    #       half of desired DLP resolution)
profile_range = [65000., 150000.]   # Lower and upper radial bounds to the
                                    #       output diffraction profile

### DiffractionCorrection
res_km = 1.0                        # Reconstruction resolution
resolution_factor = 0.75            # Factor to be multiplied to resolution
                                    #       to follow MTR86 definition
inversion_range = profile_range     # Lower and upper radial bounds to inversion
psitype = 'Fresnel4'                # Psi type
wtype = 'kbmd20'                    # Window type
sigma = 2.e-13                      # Allan deviation
fwd = False                         # Compute forward model
norm = True                         # Normalize reconstructed complex
                                    #       transmittance by window width
bfac = False                        # Use input sigma in window calculation

# ***************************

# Create new error file
err_file = sys.argv[0].split('.')[0] + time.strftime("_%Y%m%d-%H%M%S") + '.err'
fail_file = open('../output/' + err_file, 'w')
files = ['../data/'+line.strip('\n') for line in open(
         rsr_file_list,'r').readlines()]


nfiles = len(files)
init_time = time.time()

for ind in range(nfiles):

    print('\nn='+str(ind))
    rsr_file = files[ind]

    # exclude 16 kHz files
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
            planet,
            spacecraft,
            kernels,
            verbose = verbose,
            write_file = write_file
        )

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(
            rsr_inst,
            geo_inst,
            verbose = verbose,
            write_file = write_file,
            pnf_order = pnf_order,
            interact=interact
        )

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
                verbose = verbose
            )
        )

        # Invert profile for full occultation
        if dlp_inst_ing is not None:
            tau_inst = rss.DiffractionCorrection(
                dlp_inst_ing,
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

            tau_inst.history = write_history_dict(tau_inst.input_vars,
                                                  tau_inst.input_kwds,
                                                  __file__)

            if write_file:
                tau_inst.outfiles = write_output_files(tau_inst)

                rss.tools.plot_summary_doc_v2(
                    geo_inst, cal_inst, dlp_inst_ing, tau_inst
                )

        if dlp_inst_egr is not None:
            tau_inst = rss.DiffractionCorrection(
                dlp_inst_egr,
                res_km,
                rng = inversion_range,
                resolution_factor = resolution_factor,
                psitype = psitype,
                wtype = wtype,
                use_norm = norm,
                bfac = bfac,
                verbose = verbose
            )

            tau_inst.history = write_history_dict(tau_inst.input_vars,
                                                  tau_inst.input_kwds,
                                                  __file__)

            if write_file:
                tau_inst.outfiles = write_output_files(tau_inst)

                rss.tools.plot_summary_doc_v2(
                    geo_inst, cal_inst, dlp_inst_egr, tau_inst
                )

        et = time.time()
        run_time = str((et-st)/60.)
        print('File processing time (min): ' + str(run_time))
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
print('Total processing time (hrs): ' + str(total_batch_time))
fail_file.close()
