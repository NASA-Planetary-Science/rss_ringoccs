import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import matplotlib.pyplot as plt

# ***** Begin user input *****
rsr_file = "/Volumes/rmaguire002/data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX43RD.2A1"
planet = 'Saturn'
spacecraft = 'Cassini'
kernels = '../tables/Sa-TC17-V001.ker'
dr_km_desired = 0.25
rng = [65000., 150000.]
res_km = 0.75
verbose = True
write_file = False
interact = False
# ***** End user input *****

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose, write_file=write_file)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,verbose=verbose,
        write_file=write_file,fof_order=9,pnf_order=3, interact=interact)

# Create instances of DLP class with ingress (ing) and egress (egr) profiles
dlp_inst_ing, dlp_inst_egr = (
        rss.calibration.DiffractionLimitedProfile.create_dlps(
            rsr_inst, geo_inst, cal_inst, dr_km_desired, 
            write_file=write_file,verbose=verbose))

if dlp_inst_ing is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, res_km,
            rng=rng,write_file=write_file,verbose=verbose,psitype='Fresnel4')
    rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_ing, 
                                    tau_inst)
if dlp_inst_egr is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, res_km,
            rng=rng,write_file=write_file,verbose=verbose,psitype='Fresnel4')

    rss.tools.plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst_egr, tau_inst)
