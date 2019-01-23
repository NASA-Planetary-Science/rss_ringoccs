"""
e2e_run.py

Purpose:
    Example rss_ringoccs package script to produce a 1-km resolution
    profile of the Huygens ringlet for Rev007 E X43.

Notes:
    #. This program assumes that all relevant files (rsr_file and
        kernels) have already been downloaded. See p6-7 of
        rss_ringoccs: User's Guide for details.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import time

sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

# ***** Begin user input *****
rsr_file = '../data/co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX43RD.2A1'

kernels = '../examples/Rev007_meta_kernel.ker'

dr_km_desired = 0.05
res_km = 1.0

planet = 'Saturn'
spacecraft = 'Cassini'

fof_order = 9
pnf_order = 3
psitype = 'Fresnel4'

verbose = True
write_file = True
interact = False

feature_km = 117830.
feature_name = 'Huygens Ringlet'

inversion_range = [feature_km-40., feature_km+40.]


# ***** End user input *****

start_time = time.time()

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose, write_file=write_file)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst, verbose=verbose,
        write_file=write_file, fof_order=fof_order, pnf_order=pnf_order,
        interact=interact)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst_ing, dlp_inst_egr = (
        rss.calibration.DiffractionLimitedProfile.create_dlps(
            rsr_inst, geo_inst, cal_inst, dr_km_desired, write_file=write_file,
            verbose=verbose))

# Invert profile for full occultation
if dlp_inst_ing is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, res_km,
            rng=inversion_range, write_file=write_file, verbose=verbose,
            psitype=psitype)
if dlp_inst_egr is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, res_km,
            rng=inversion_range, write_file=write_file, verbose=verbose,
            psitype=psitype)


# Plot Huygens ringlet
#   first row of three: uncorrected optical depth, power, and phase
#   second row of three: corrected optical depth, power, and phase
fig, axes = plt.subplots(2,3, figsize=(8,5))

rho_km = tau_inst.rho_km_vals - feature_km
recon_power = tau_inst.power_vals
raw_power = tau_inst.p_norm_vals

raw_tau = -tau_inst.mu_vals * np.log(tau_inst.p_norm_vals)
recon_tau = tau_inst.tau_vals

raw_phase = tau_inst.phase_rad_vals
recon_phase = tau_inst.phase_vals
xtitle = '$\\rho$ - ' + str(feature_km) + ' (km)'

fig.suptitle('Rev007E X43 ' + feature_name + ' at ' + str(res_km) + ' km resolution', fontsize=14)

axes[0,0].plot(rho_km, raw_power)
axes[0,0].set_ylabel('Raw Power')
axes[0,0].set_xlabel(xtitle)

axes[0,1].plot(rho_km, raw_tau)
axes[0,1].set_ylabel('Raw Optical Depth')
axes[0,1].set_xlabel(xtitle)

axes[0,2].plot(rho_km, raw_phase)
axes[0,2].set_ylabel('Raw Phase')
axes[0,2].set_xlabel(xtitle)

axes[1,0].plot(rho_km, recon_power)
axes[1,0].set_ylabel('Reconstructed Power')
axes[1,0].set_xlabel(xtitle)

axes[1,1].plot(rho_km, recon_tau)
axes[1,1].set_ylabel('Reconstructed Optical Depth')
axes[1,1].set_xlabel(xtitle)

axes[1,2].plot(rho_km, recon_phase)
axes[1,2].set_ylabel('Reconstructed Phase')
axes[1,2].set_xlabel(xtitle)

end_time = time.time()
plt.tight_layout()
plt.subplots_adjust(top=0.93)

print('Total run time (min): ', str((end_time-start_time)/60.))

plt.show()

"""
Revisions:
"""
