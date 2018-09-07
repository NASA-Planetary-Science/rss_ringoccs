"""

Rev007EX43_HuygensRinglet_test.py

Purpose: rss_ringoccs package sample run to produce a 1-km resolution profile
         of the Huygens ringlet using Rev007E X43 data. To speed up the process,
         pre-computed files will be used (f_resid_fit_parameters.p,
         power_norm_fit_parameters.p, freq_offset_file.txt). These files can
         be found in rss_ringoccs/output/rev7E_X43_e2e_output/.

Notes:
    [1] This program assumes that all relevant files (rsr_file and
        kernels) have already been downloaded. See p XX of
        rss_ringoccs: User's Guide for details.
    [2] This program does not produce any output files.

Revisions:
    2018 Aug 21 - jfong - original
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
#rsr_file = '../data/co-s-rss-1-sroc1-v10/cors_0108/sroc1_123/rsr/s10sroi2005123_0230nnnx43rd.2a2'
#rsr_file = '../../../../data/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2'
rsr_file = '../../../../data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
kernels_list_file = '../tables/Rev007_list_of_kernels.txt'
kernels = '../tables/Rev007_meta_kernel.ker'
kernels_dir = '../kernels/'

output_directory = '../output/rev7E_X43_e2e_output/'
#output_directory = '../../../../../jfong/rss_ringoccs/output/RSS_loop_v2_output/Rev007/E/Rev007E_RSS_2005_123_X43_E/'
freq_offset_file = output_directory + 'freq_offset_file.txt'
f_resid_fit_parameters_file = output_directory + 'f_resid_fit_parameters.p'
power_norm_fit_parameters_file = (output_directory
    + 'power_norm_fit_parameters.p')

outfig = '../../../../../jfong/rss_ringoccs/output/figs/users_guide_Huygens_20180822.ps'



f_USO = 8427222034.34050
dr_km_desired = 0.25
res_km = 1.0

verbose = True
planet = 'Saturn'
spacecraft = 'Cassini'

Huygens_xrange = [117800., 117860.] 

# ***** End user input *****

start_time = time.time()
def read_f_resid_fit_parameters(f_resid_fit_parameters_file):

    file_object = open(f_resid_fit_parameters_file, 'rb')
    fit_param_dict = pickle.load(file_object)
    k = fit_param_dict['k']
    spm_include = fit_param_dict['spm_include']
    return k, spm_include

def read_power_norm_fit_parameters(power_norm_fit_parameters_file):

    file_object = open(power_norm_fit_parameters_file, 'rb')
    fit_param_dict = pickle.load(file_object)
    k = fit_param_dict['k']
    freespace_spm = fit_param_dict['freespace_spm']
    knots_spm = fit_param_dict['knots_spm']
    return k, freespace_spm, knots_spm

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose)

# Read in frequency offset data from pre-computed file
freq_offset_file_vals = np.loadtxt(freq_offset_file)
f_spm = freq_offset_file_vals[:, 0]
f_offset = freq_offset_file_vals[:, 1]

# Read in residual frequency fit parameters from pre-computed file
k_f_resid, spm_include = read_f_resid_fit_parameters(
        f_resid_fit_parameters_file)

fit_inst = rss.calibration.FreqOffsetFit(rsr_inst, geo_inst, f_spm,
        f_offset, f_USO, poly_order=k_f_resid, spm_include=spm_include,
        USE_GUI=False, verbose=verbose)

# Get corrected I's and Q's
spm_vals, IQ_c = fit_inst.get_IQ_c()

# Create instance with normalized power
norm_inst = rss.calibration.Normalization(spm_vals, IQ_c, geo_inst, rsr_inst,
    verbose=verbose)

# Read in power normalization fit parameters from pre-computed fole
k_power_norm, freespace_spm, knots_spm = read_power_norm_fit_parameters(
        power_norm_fit_parameters_file)
spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
        freespace_spm=freespace_spm, knots_spm=knots_spm,
        spline_order=k_power_norm, USE_GUI=False, verbose=verbose)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(fit_inst, norm_inst, geo_inst,
        verbose=verbose)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst = rss.calibration.NormDiff(rsr_inst, dr_km_desired, geo_inst,
        cal_inst, verbose=verbose)

# Invert profile for only Huygens ringlet
tau_inst = rss.diffcorr.DiffractionCorrection(dlp_inst, res_km,
        rng=Huygens_xrange, verbose=verbose)

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
plt.savefig(outfig)

print('Total run time: ', end_time-start_time)
plt.show()

