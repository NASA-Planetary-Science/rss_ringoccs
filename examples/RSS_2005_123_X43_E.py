"""

RSS_2005_123_X43_E.py

Purpose: Example of a possible end-to-end run. Produces a calibration file and
         obs file for rev7E X43

NOTE: This was made for the 1kHz rev7E X43 file, which is unfortunately not
      available online. To run this on the 16kHz file, go to the cors_0105
      directory online, and download "s10sroe2005123_0740nnnx43rd.2a2". Then
      use that file as "rsr_file", and you probably want to specify the
      "decimate_16khz_to_1khz" keyword in "rsr_inst.get_IQ" and the norm_inst
      definition, because otherwise it will take a really long time to run.

Revisions:
      example_end_to_end_and_quick_look_script.py
   Mar 20 2018 - gsteranka - Original version
      RSS_2005_123_X43_E.py
   Apr 03 2018 - gsteranka - Copied from original to make it clear that
                             this produces Rev7E X43 stuff for input and
                             output
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
rsr_file = '../../../data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
geo_file = '../input/Rev7E/RSS_2005_123_X43_E_GEO.TAB'
cal_file = '../input/Rev7E/RSS_2005_123_X43_E_CAL.TAB'
obs_file = '../input/Rev7E/RSS_2005_123_X43_E_OBS.TAB'
out_file = '../output/Rev7E/RSS_2005_123_X43_E_OUT_01KM.TAB'

# Change to USO frequency on date closest to occultation
f_uso = 8427222034.34050

# Kernels you need and where to them
kernels_dir = '../../../kernels~/'
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
res_km = 1.0
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
