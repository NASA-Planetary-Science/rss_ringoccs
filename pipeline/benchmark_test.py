"""

Rev007EX43_time_test.py

Purpose: rss_ringoccs package sample run to produce a 1-km resolution profile
         of the full Rev007E X43 occultation. Since this is intended as a time
         test, no pre-computed files are used.

Notes:
    [1] This program assumes that all relevant files (rsr_file and
        kernels) have already been downloaded. See p XX of
        rss_ringoccs: User's Guide for details.
    [2] This program does not produce any output files.

Revisions:
        Rev007EX43_HuygensRinglet_test.py
    2018 Aug 21 - jfong - original
        Rev007EX43_time_test.py
    2018 Aug 31 - gsteranka - Edited to not use any pre-computed files, use
                              default fits for both cases rather than GUI, and
                              to not plot anything
"""

import numpy as np
import pdb
import sys
import time

sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

# ***** Begin user input *****
rsr_file = ('../data/co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/'
    + 's10sroe2005123_0740nnnx43rd.2a2')
kernels = '../tables/Rev007_meta_kernel.ker'

f_USO = 8427222034.34050
dr_km_desired = 0.05
res_km = 0.100

verbose = True
planet = 'Saturn'
spacecraft = 'Cassini'

inversion_range = [70000, 145000]

# ***** End user input *****

start_time = time.time()

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose)

f_spm, f_offset, f_offset_history = rss.calibration.calc_freq_offset(rsr_inst)

# Use default residual frequency fit
fit_inst = rss.calibration.FreqOffsetFit(rsr_inst, geo_inst, f_spm,
        f_offset, f_USO, USE_GUI=False, verbose=verbose)

# Get corrected I's and Q's
spm_vals, IQ_c = fit_inst.get_IQ_c()

# Create instance with normalized power
norm_inst = rss.calibration.Normalization(spm_vals, IQ_c, geo_inst, rsr_inst,
    verbose=verbose)

# Use default power fit
spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(USE_GUI=False,
    verbose=verbose)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(fit_inst, norm_inst, geo_inst,
        verbose=verbose)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst = rss.calibration.NormDiff(rsr_inst, dr_km_desired, geo_inst,
        cal_inst, verbose=verbose)

# Invert profile for full occultation
tau_inst = rss.diffcorr.DiffractionCorrection(dlp_inst, res_km,
        rng=inversion_range, verbose=verbose)

end_time = time.time()
print('Total run time: ', end_time-start_time)

pdb.set_trace()
