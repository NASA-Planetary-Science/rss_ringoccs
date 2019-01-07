import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import pdb
import time

# ***** Begin user input *****
rsr_file_list = [
        "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNK34RD.1B1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNS43RD.2B1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX34RD.1A1",
        "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX43RD.2A1"
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0228NNNX14RD.2A1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0230NNNK34RD.1B1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0230NNNS14LD.1N1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0230NNNS14RD.2B1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0230NNNS43RD.2B1",
        # "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROI2005123_0230NNNX34RD.1A1"
]

kernels = '../tables/Sa-TC17-V001.ker'
planet = 'Saturn'
spacecraft = 'Cassini'
dr_km_desired = 0.25
res_km = 0.75
inversion_range = [65000., 150000.]
psitype='Fresnel'

#USE_GUI = True
write_file = True
verbose = True
# ***** End user input *****

start_time = time.time()
for rsr_file in rsr_file_list:
    # Create instance with rsr file contents
    rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)
    
    # Create instance with geometry parameters
    geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
            verbose=verbose, write_file=write_file)
    
    # Create instance with calibrated data
    cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
            verbose=verbose,
            write_file=write_file)
    
    
    
    # Create instance with diffraction-limited profile and other
    #   inputs needed for diffraction correction
    dlp_inst_ing, dlp_inst_egr = rss.calibration.DiffractionLimitedProfile.create_dlps(rsr_inst, geo_inst, cal_inst, dr_km_desired, verbose=verbose, write_file=write_file)

    if dlp_inst_ing:
        tau_inst_ing = rss.diffrec.DiffractionCorrection(dlp_inst_ing,
                res_km, rng=inversion_range, verbose=verbose, write_file=write_file,
                psitype=psitype)
    
        
    if dlp_inst_egr:
        tau_inst_egr = rss.diffrec.DiffractionCorrection(dlp_inst_egr,
                res_km, rng=inversion_range, verbose=verbose, write_file=write_file,
                psitype=psitype)
