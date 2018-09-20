"""

benchmark_run.py

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
    2018 Sep 20 - jfong - updated to work with modified e2e structure
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

dr_km_desired = 0.05
res_km = 0.10

file_search = False
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

end_time = time.time()
print('Total run time: ', end_time-start_time)

pdb.set_trace()
