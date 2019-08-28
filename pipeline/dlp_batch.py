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

    if rsr_file[-1] == '5' :
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

        # Create spectrogram of scattered signal
        #rss.scatter.Scatter(rsr_inst,geo_inst,cal_inst,stack=True,nstack=int(9))

    # allow KeyboardInterrupt
    except KeyboardInterrupt:
        sys.exit()
    # output errors to file
    except:
        tb = traceback.format_exc()
        print(tb)
        fail_file.write('-'*72+'\n')
        fail_file.write('  '+dir+'\n'+'-'*72+'\n\n'+ 'n='+
                         str(n)+'\n'+ tb+'\n')
    en = time.time()

    print('Initial processing time:   ',round(en-st,3))
