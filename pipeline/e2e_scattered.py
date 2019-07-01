import numpy as np
import sys
import time
import traceback
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import e2e_batch_args as args
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

rho_limits = [6.5e4,1.4e5]

### Read in file names and Maxwell values
files = open('../tables/rsr_1kHz_files_before_USO_failure.txt','r').readlines()

# do for all files
for file in files:
    st = time.time()
    # RSR file name
    rsr_file = '../data/'+file.strip('\n')
    print(rsr_file)

    # exclude 16 kHz files
    if rsr_file[-1] == '2' and not args.with16:
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        continue

    try:
        # Create instance with rsr file contents
        rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=args.verbose,
                decimate_16khz_to_1khz=args.decimate_16khz_to_1khz)

        # Create instance with geometry parameters
        geo_inst = rss.occgeo.Geometry(rsr_inst, args.planet, args.spacecraft,
                args.kernels, verbose=args.verbose, write_file=args.write_file)

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
                verbose=args.verbose, write_file=args.write_file,
                pnf_order=args.pnf_order, interact=args.interact)

        # Create instance with diffraction-limited profile and other
        #   inputs needed for diffraction correction
        dlp_inst_ing, dlp_inst_egr = (
                rss.calibration.DiffractionLimitedProfile.create_dlps(
                    rsr_inst, geo_inst, cal_inst, args.dr_km_desired,
                    profile_range = args.profile_range,
                    write_file=args.write_file, verbose=args.verbose))
        # Invert profile for full occultation
        if dlp_inst_ing is not None:
            tau_inst = (rss.diffrec.DiffractionCorrection(
                    dlp_inst_ing, args.res_km,
                    rng=args.inversion_range, res_factor=args.res_factor,
                    psitype=args.psitype, wtype=args.wtype, fwd=args.fwd,
                    norm=args.norm, bfac=args.bfac, write_file=args.write_file,
                    verbose=args.verbose))
            rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_ing,
                    tau_inst)
        if dlp_inst_egr is not None:
            tau_inst = (rss.diffrec.DiffractionCorrection(
                    dlp_inst_egr, args.res_km,
                    rng=args.inversion_range, res_factor=args.res_factor,
                    psitype=args.psitype, wtype=args.wtype, fwd=args.fwd,
                    norm=args.norm, bfac=args.bfac, write_file=args.write_file,
                    verbose=args.verbose))
            rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_egr,
                    tau_inst)

        # Create spectrogram of scattered signal
        rss.scatter.Scatter(rsr_inst,geo_inst,cal_inst,stack=True,nstack=int(9))
    except KeyboardInterrupt:
        sys.exit()
    except:
        tb = traceback.format_exc()
        print(tb)
    en = time.time()
    print(en-st)
