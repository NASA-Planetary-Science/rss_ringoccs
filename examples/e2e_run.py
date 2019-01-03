"""
Purpose: Example 'End-to-End' script; refer to pXX of User's Guide.
"""
import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

###
### ***** Begin user input *****
###
rsr_file = ('../data/co-s-rss-1-sroc11-v10/cors_0753/SROC11_026/RSR/'
            + 'S57SROE2010027_0005NNNK55RD.1B1')
kernels = '../tables/Sa-TC17-V001.ker'
planet = 'Saturn'
spacecraft = 'Cassini'
# DLP resolution
dr_km_desired = 0.05
# reconstruction resolution
res_km = 0.75
# reconstruction radius range
rng = [6e4,1.5e5]
# output keyword options
write_file = True
verbose = False

###
### ***** End user input *****
###

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose, decimate_16khz_to_1khz=False)
# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose, write_file=write_file)
# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,verbose=verbose,
        write_file=write_file,fof_order=9,pnf_order=3)
# Create instances of DLP class with ingress (ing) and egress (egr) profiles
dlp_inst_ing, dlp_inst_egr = rss.calibration.DiffractionLimitedProfile.create_dlps(rsr_inst, geo_inst,
                                            cal_inst, dr_km_desired,write_file=write_file,verbose=verbose)
# Create instance of diffraction reconstruction for each DLP instance
if dlp_inst_ing is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, res_km,
            rng=rng,write_file=write_file,verbose=verbose,psitype='Fresnel4')
if dlp_inst_egr is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, res_km,
            rng=rng,write_file=write_file,verbose=verbose,psitype='Fresnel4')
