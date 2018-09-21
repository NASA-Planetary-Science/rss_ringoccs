"""

e2e_run_with_files.py

Purpose: Example rss_ringoccs package script to produce a 500-m resolution
         profile of the Huygens ringlet using pre-computed intermediate files
         for Rev007 E X43.

Notes:
    [1] This program assumes that all relevant files (rsr_file and
        kernels) have already been downloaded. See p6-7 of
        rss_ringoccs: User's Guide for details.

Revisions:
    2018 Sep 19 - jfong - original
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
import time

sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

# ***** Begin user input *****
#rsr_file = '../data/co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2'
rsr_file = '/Volumes/jfong001/Research/TC2017/data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
kernels_list_file = 'Rev007_list_of_kernels.txt'
kernels = 'Rev007_meta_kernel.ker'
kernels_dir = '../kernels/'

dr_km_desired = 0.25
res_km = 0.5

file_search = True
verbose = True
planet = 'Saturn'
spacecraft = 'Cassini'

inversion_range = 'all'

# ***** End user input *****

start_time = time.time()

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
        file_search=file_search, verbose=verbose)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst = rss.calibration.NormDiff(rsr_inst, geo_inst, cal_inst,
        dr_km_desired, verbose=verbose)

# Invert profile for full occultation
tau_inst = rss.diffcorr.DiffractionCorrection(dlp_inst, res_km,
        rng=inversion_range, verbose=verbose)

# Plot Huygens ringlet
#   first row of three: uncorrected optical depth, power, and phase
#   second row of three: corrected optical depth, power, and phase
fig, axes = plt.subplots(2,3, figsize=(8,5))

rho_km = tau_inst.rho_km_vals - 117830.
recon_power = tau_inst.power_vals
raw_power = tau_inst.p_norm_vals

raw_tau = -tau_inst.mu_vals * np.log(tau_inst.p_norm_vals)
recon_tau = tau_inst.tau_vals

raw_phase = tau_inst.phase_rad_vals
recon_phase = tau_inst.phase_vals
xtitle = '$\\rho$ - 117830 (km)'

fig.suptitle('Rev007E X43 Huygens Ringlet at 1km resolution', fontsize=14)

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

print('Total run time: ', end_time-start_time)


