"""

example_end_to_end_and_quick_look_script.py

Purpose: Example of a possible end-to-end run. Produces a calibration file

Revisions:
      example_end_to_end_and_quick_look_script.py
   Mar 20 2018 - gsteranka - Original version
"""

# Use sys.path.append here to add directory where rest of code is.
# Typically "sys.path.append('../src/')
import sys
sys.path.append('../src')

import numpy as np

from rsr_reader import RSRReader
from freq_offset import calc_freq_offset
from calc_geometry_v2 import Geometry
from freq_offset_fit import FreqOffsetFit
from power_normalization import Normalization
from make_cal_file import make_cal_file
from make_geo_file import make_geo_file
from norm_diff_class import NormDiff
from diffraction_correction import diffraction_correction

# Change to location of rsr file, and where to write out geo
# and cal files
rsr_file = '/your/path/name/here/example_rsr_file'
geo_file = '/your/path/name/here/example_geo_file.tab'
cal_file = '/your/path/name/here/example_cal_file.tab'
obs_file = '/your/path/name/here/example_obs_file.tab'
out_file = '/your/path/name/here/example_out_file.tab'

# Change to USO frequency on date closest to occultation
f_uso = 8427222034.34050

# Kernels you need and where to them
kernels_dir = '/your/path/name/here/'
kernel_files = ['de421.bsp',
                '050606R_SCPSE_05114_05132.bsp',
                'earthstns_itrf93_040916.bsp',
                'earth_720101_070426.bpc',
                'cpck26Feb2009.tpc',
                'naif0012.tls']
kernels = [kernels_dir + this_kernel for this_kernel in kernel_files]

# NormDiff input in quick-look steps
dr_km_desired = 0.25
rho_km_range = [70000, 155000]
res_km = 'I do nothing yet'
window_type = 'I do nothing yet'

# **END OF USER INPUT**

rsr_inst = RSRReader(rsr_file)
(spm_vals, IQ_m) = rsr_inst.get_IQ()

(f_spm_offset, f_offset) = calc_freq_offset(spm_vals, IQ_m)

geo_inst = Geometry(rsr_inst, kernels)

fit_inst = FreqOffsetFit(rsr_inst, geo_inst, f_spm_offset, f_offset,
                         f_uso, kernels)
(dummy_spm, IQ_c) = fit_inst.get_IQ_c()

norm_inst = Normalization(spm_vals, IQ_c, geo_inst)
spm_fit = np.arange(spm_vals[0], spm_vals[-1], 1.0)
spline_fit = norm_inst.get_spline_fit(spm_fit)

make_geo_file(geo_inst, geo_file)
make_cal_file(cal_file, fit_inst, norm_inst, geo_inst)

# Beginning of quick-look steps - from here on out, you can run these steps
# without running the above, since the geo and cal files take care of the
# above steps
norm_diff_inst = NormDiff(rsr_file, dr_km_desired, rho_km_range, res_km,
                          window_type, geo_file, cal_file)
norm_diff_inst.save_obs(obs_file)

# After making obs file, you can repeat this step as much as needed
# without repeating any of the above
diffraction_correction(obs_file, res_km, out_file)
