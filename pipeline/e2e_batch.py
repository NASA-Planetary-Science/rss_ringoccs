import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import traceback
import time
###
### ***** BEGIN USER INPUT *****
###
# specify the directory of the RSR data files
mpath = '../data/'
# specify the meta-kernel files
kernels = '../tables/e2e_kernels.ker'
# specify the planet
planet = 'Saturn'
# specify the spacecraft
spacecraft = 'Cassini'
# specify whether to process 16 kHz files
# WARNING: THESE REQUIRE A MINIMUM OF
# 2.9 GHZ CPU, 16 GB RAM TO PROCESS
with16 = False
# specify the DLP radial sampling rate in km
dr_km_desired = 0.05
# specify the TAU processing resolution in km
res_km = 0.5
# specify the radial range over which to
# reconstruct the optical depth profile
inversion_range = [6e4,1.5e5]
# specify the frequency offset fit order
fof_order = 9
# specify the power normalization order
pnf_order = 3
# specify the method of phase approximation
psitype='Fresnel4'
# specify whether to output results files
write_file = True
# specify whether to output progress to terminal
verbose = False

err_file = sys.argv[0].split('.')[0] + '.err'
fail_file = open('../output/' + err_file, 'w')
#
###
### ***** END USER INPUT *****
###

# import files
files = [mpath+line.strip('\n') for line in open('../tables/'+
    'list_of_rsr_files_before_USO_failure_to_dl_v2.txt','r').readlines()]

# files with bad headers (to exclude)
skips = [mpath+'co-s-rss-1-sroc8-v10/cors_0745/SROC8_239/RSR/S43SROI2008239_1410NNNX63RD.1A1',
mpath+'co-s-rss-1-sroc8-v10/cors_0745/SROC8_239/RSR/S43SROI2008239_1410NNNS63RD.1B1']

init_time = time.time()

for rsr_file in files:

    # exclude 16 kHz files
    if rsr_file[-1] == '2' and not with16:
        print('SKIPPING 16 KHZ FILE: ' + rsr_file)
        continue
    # exclude files with bad headers
    if rsr_file in skips:
        print('SKIPPING FILE WITH BAD HEADER: '+rsr_file)
        continue
    try:
        st = time.time()
    
        # print RSR file
        print(rsr_file)
    
        # Create instance with rsr file contents
        rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose,
                decimate_16khz_to_1khz=with16)
    
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
