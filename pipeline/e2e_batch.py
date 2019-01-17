import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

import traceback

import time

# ***** Begin user input *****

mpath = '../data/'

kernels = '../tables/e2e_kernels.ker'
planet = 'Saturn'
spacecraft = 'Cassini'
dr_km_desired = 0.25
res_km = 1.0
inversion_range = [6e4,1.5e5]
fof_order = 9
pnf_order = 3

fail_file = open('final_failures.txt', 'w')

write_file = True
verbose = False
# ***** End user input *****

files = [mpath+line.strip('\n') for line in open('../tables/list_of_rsr_files_before_USO_failure_to_dl_v2.txt','r').readlines()]
init_time = time.time()
nfiles = len(files)
for ind in range(nfiles):

    print('n=' + str(ind))
    rsr_file = files[ind]

    if rsr_file[-1] == '2':
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        continue

    try:
        st = time.time()

        # make RSR file
        print(rsr_file)

        # Create instance with rsr file contents
        rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose,
                decimate_16khz_to_1khz=False)

        # Create instance with geometry parameters
        geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
                verbose=verbose, write_file=write_file)

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
                verbose=verbose, write_file=write_file,
                fof_order=fof_order, pnf_order=pnf_order)

        # Create instaces of DLP class for one or both directions
        # depending on occultation type (diametric or chord, respectively)
        dlp_inst_ing, dlp_inst_egr = (
                rss.calibration.DiffractionLimitedProfile.create_dlps(
                    rsr_inst, geo_inst, cal_inst, dr_km_desired,
                    write_file=write_file, verbose=verbose))

        # Create instances of the diffraction-reconstructed profile
        # for existing instances of the DLP
        if dlp_inst_ing is not None:
            tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, res_km,
                    write_file=write_file, rng=inversion_range,
                    verbose=verbose, psitype=psitype)
            rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_ing,
                    tau_inst)
        if dlp_inst_egr is not None:
            tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, res_km,
                    write_file=write_file, rng=inversion_range,
                    verbose=verbose, psitype=psitype)
            rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_egr,
                    tau_inst)


        # if verbose, let user know this file processed successfuly
        et = time.time()
        run_time = str((et-st)/60.)
        print('File processing time (min): ' + str(run_time))
        #if verbose:
        #    print('-'*48+'\n')
        #    print('  FILE '+file+' PROCESSED\n'+'-'*48+'\n')


    except:
        # if verbose, let user know this file failed to process and why
        #if verbose:
        tb = traceback.format_exc()
        fail_file.write('-'*96+'\n')
        print(tb)
        fail_file.write('  '+rsr_file+'\n'+'-'*96+'\n\n'+ 'n='+
                         str(ind)+'\n'+ tb+'\n')

fail_file.close()

final_time = time.time()
total_batch_time = (final_time - init_time)/60./60.
print('Total processing time (hrs): ' + str(total_batch_time))
