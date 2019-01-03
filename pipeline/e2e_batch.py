import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

# ***** Begin user input *****
mpath = '/Volumes/sflury001/Research/TC2017/sflury/rss_ringoccs/data'

kernels = '../tables/Sa-TC17-V001.ker'
planet = 'Saturn'
spacecraft = 'Cassini'
dr_km_desired = 0.05
res_km = 0.75
inversion_range = [6e4,1.5e5]

write_file = True
verbose = False
# ***** End user input *****

files = [line.strip('\n') for line in open('../tables/list_of_rsr_files_before_USO_failure_to_dl_v2.txt','r').readlines()]

for file in files:

    try:

        # make RSR file
        rsr_file = mpath+'/'+file
        print('\n'+rsr_file+'\n')

        # Create instance with rsr file contents
        rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose, decimate_16khz_to_1khz=False)

        # Create instance with geometry parameters
        geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
                verbose=verbose, write_file=write_file)

        # Create instance with calibrated data
        cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,verbose=verbose,
                write_file=write_file,fof_order=9,pnf_order=3)

        # Create instaces of DLP class for one or both directions
        # depending on occultation type (diametric or chord, respectively)
        dlp_inst_ing, dlp_inst_egr = rss.calibration.DiffractionLimitedProfile.create_dlps(rsr_inst, geo_inst,
                                                    cal_inst, dr_km_desired,write_file=write_file,verbose=verbose)

        # Create instances of the diffraction-reconstructed profile
        # for existing instances of the DLP
        if dlp_inst_ing is not None:
            tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, res_km,
                    write_file=write_file, rng=[6e4,1.5e5], verbose=verbose,psitype='Fresnel4')
        if dlp_inst_egr is not None:
            tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, res_km,
                    write_file=write_file, rng=[6e4,1.5e5], verbose=verbose,psitype='Fresnel4')

        # if verbose, let user know this file processed successfuly
        if verbose:
            print('-'*48+'\n')
            print('  FILE '+file+' PROCESSED\n'+'-'*48+'\n')

    except:
        # if verbose, let user know this file failed to process and why
        if verbose:
            tb = traceback.format_exc()
            print('-'*48+'\n')
            print('  FILE '+file+' FAILED\n'+'-'*48+'\n\n'+tb+'\n')
