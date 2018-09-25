"""

e2e_run.py

Purpose: Sample 'End-to-End' script; refer to pXX of User's Guide.

Revisions:
    2018 Sep 17 - jfong - original
"""

import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import pdb
import time

# ***** Begin user input *****
rsr_file = '../data/co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2'
kernels = 'Rev007_meta_kernel.ker'
planet = 'Saturn'
spacecraft = 'Cassini'
dr_km_desired = 0.05
res_km = 0.10
inversion_range = 'all'

USE_GUI = True
file_search = False
write_file = True
verbose = True
# ***** End user input *****

start_time = time.time()

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

# Create instance with geometry parameters
geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft, kernels,
        verbose=verbose, write_file=write_file)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst,
        file_search=file_search, USE_GUI=USE_GUI, verbose=verbose,
        write_file=write_file)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst = rss.calibration.NormDiff(rsr_inst, geo_inst,
        cal_inst, dr_km_desired, verbose=verbose, write_file=write_file)

# Invert profile for full occultation
tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst, res_km,
        rng=inversion_range, verbose=verbose, write_file=write_file)
end_time = time.time()
print('Computation time: ' + str(end_time-start_time))
pdb.set_trace()
