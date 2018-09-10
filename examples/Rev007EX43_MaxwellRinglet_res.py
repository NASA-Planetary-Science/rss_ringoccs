"""

Rev007EX43_MaxwellRinglet_res.py

Purpose: rss_ringoccs package sample run to produce a 750m, 500m, 250m, 100m
         resolution profile of the Maxwell ringlet using Rev007E X43 data.
         To speed up the process, pre-computed files will be used 
         (f_resid_fit_parameters.p, power_norm_fit_parameters.p, 
         freq_offset_file.txt). These files can be found in 
         rss_ringoccs/output/rev7E_X43_e2e_output/.

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
#rsr_file = '../data/co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2'
rsr_file = '/Volumes/jfong001/Research/TC2017/data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
kernels_list_file = '../tables/Rev007_list_of_kernels.txt'
kernels = '../tables/Rev007_meta_kernel.ker'
kernels_dir = '../kernels/'

output_directory = '../output/rev7E_X43_e2e_output/'
freq_offset_file = output_directory + 'freq_offset_file.txt'
f_resid_fit_parameters_file = output_directory + 'f_resid_fit_parameters.p'
power_norm_fit_parameters_file = (output_directory
    + 'power_norm_fit_parameters.p')

f_USO = 8427222034.34050
dr_km_desired = 0.25

verbose = True
planet = 'Saturn'
spacecraft = 'Cassini'

Maxwell_xrange = [87400., 87600.] 
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

# Invert profile for only Maxwell ringlet, at 750m, 500m, 250m, and 100m
print('Reconstructing at 750m...')
tau_inst_750m = rss.diffcorr.DiffractionCorrection(dlp_inst, 0.75,
        rng=Maxwell_xrange, verbose=verbose)

print('Reconstructing at 500m...')
tau_inst_500m = rss.diffcorr.DiffractionCorrection(dlp_inst, 0.50,
        rng=Maxwell_xrange, verbose=verbose)

print('Reconstructing at 250m...')
tau_inst_250m = rss.diffcorr.DiffractionCorrection(dlp_inst, 0.25,
        rng=Maxwell_xrange, verbose=verbose)

print('Reconstructing at 100m...')
tau_inst_100m = rss.diffcorr.DiffractionCorrection(dlp_inst, 0.10,
        rng=Maxwell_xrange, verbose=verbose)
# Plot Maxwell ringlet
#   first row of three: uncorrected optical depth, power, and phase
#   second row of three: corrected optical depth, power, and phase
fig, axes = plt.subplots(4,1, figsize=(8.5,11), sharex=True)
plt.subplots_adjust(hspace=0)

rho_750 = tau_inst_750m.rho_km_vals - 87515.
tau_750 = tau_inst_750m.tau_vals

rho_500 = tau_inst_500m.rho_km_vals - 87515.
tau_500 = tau_inst_500m.tau_vals

rho_250 = tau_inst_250m.rho_km_vals - 87515.
tau_250 = tau_inst_250m.tau_vals

rho_100 = tau_inst_100m.rho_km_vals - 87515.
tau_100 = tau_inst_100m.tau_vals


fig.suptitle('Rev007E X43 Maxwell Ringlet Optical Depth \nReconstruction Resolution Comparison', fontsize=13)
axes[0].plot(rho_750, tau_750, label='750m')
axes[1].plot(rho_500, tau_500, label='500m')
axes[2].plot(rho_250, tau_250, label='250m')
axes[3].plot(rho_100, tau_100, label='100m')

for ax in axes:
    ax.grid(True)
    ax.set_ylabel('$\\tau$')
    ax.set_xlabel('$\\rho$ - 87515 (km)')
#    ax.set_xlim([87470.,87560.])
    ax.set_xlim([-45.,45.])
    ax.legend()



end_time = time.time()
print('Total run time: ', end_time-start_time)
plt.show()

