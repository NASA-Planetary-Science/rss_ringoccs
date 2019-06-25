import e2e_batch_args as args
import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import traceback
import time

# Create new error file
err_file = sys.argv[0].split('.')[0] + time.strftime("_%Y%m%d-%H%M%S") + '.err'
fail_file = open('../output/' + err_file, 'w')
files = [args.mpath+line.strip('\n') for line in open(
                args.rsr_file_list,'r').readlines()]



nfiles = len(files)
init_time = time.time()

for ind in range(nfiles):
    
    print('\nn='+str(ind))
    rsr_file = files[ind]

    # exclude 16 kHz files
    if rsr_file[-1] == '2' and not args.with16:
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        continue

    try:
        st = time.time()
    
        # print RSR file
        print(rsr_file)
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
        
        et = time.time()
        run_time = str((et-st)/60.)
        print('File processing time (min): ' + str(run_time))
    except KeyboardInterrupt:
        sys.exit()
    except:
        tb = traceback.format_exc()
        fail_file.write('-'*48+'\n')
        print(tb)
        fail_file.write('  '+rsr_file+'\n'+'-'*36+'\n\n'+ 'n='+
                         str(ind)+'\n'+ tb+'\n')


final_time = time.time()
total_batch_time = (final_time - init_time)/60./60.
print('Total processing time (hrs): ' + str(total_batch_time))
fail_file.close()
