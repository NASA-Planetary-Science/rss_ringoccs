"""

example_end_to_end_script.py

Purpose: Example of a possible end-to-end run. Produces a calibration file

Revisions:
      example_end_to_end_script.py
   Mar 20 2018 - gsteranka - Original version
"""

# Use sys.path.append here to add directory where rest of code is.
# Typically "sys.path.append('../src/')
import sys
sys.path.append('../src')

import numpy as np

from rsr_reader import RSRReader
from freq_offset import calc_freq_offset
#from geo_file_into_instance import geo_file_into_instance
from calc_geometry_v2 import Geometry
from freq_offset_fit import FreqOffsetFit
from power_normalization import Normalization
from make_cal_file import make_cal_file
from make_geo_file import make_geo_file
from norm_diff_class import NormDiff

rsr_file = 'rev007E_X43/S10EAOE2005_123_0740NNNX43D.2A1'
#geo_file = 'rev007E_X43/rev007_E_X_geo.tab'
geo_file = 'rev007E_X43/test_geo.tab'

f_uso = 8427222034.34050
kernels_dir = '../../../../kernels~/'
kernel_files = ['de421.bsp',
                '050606R_SCPSE_05114_05132.bsp',
                'earthstns_itrf93_040916.bsp',
                'earth_720101_070426.bpc',
                'cpck26Feb2009.tpc',
                'naif0012.tls']
kernels = [kernels_dir + this_kernel for this_kernel in kernel_files]

#cal_file = 'rev007E_X43/rev007_E_X43_cal_example.tab'
cal_file = 'rev007E_X43/test_cal.tab'

dr_km_desired = 0.25
rho_km_range = [70000, 155000]
res_km = 'I do nothing yet'
window_type = 'I do nothing yet'

obs_file = 'rev007E_X43/test_obs.tab'

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

# Put diffraction correction code here
