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

import numpy as np
import pdb
import sys

# Use sys.path.append here to add directory where rest of code is.
#     Typically "sys.path.append('../src/')"
sys.path.append('../src')
from rsr_reader import RSRReader
from calc_freq_offset import calc_freq_offset
from calc_occ_geometry import Geometry
from freq_offset_fit import FreqOffsetFit
from power_normalization import Normalization
from make_cal_file import make_cal_file
from make_geo_file import make_geo_file
from norm_diff_class import NormDiff
from make_dlp_file import make_dlp_file
from make_cal_inst import MakeCalInst
from calibration_class import Calibration

# Change to location of rsr file, and where to write out geo
# and cal files
rsr_file = '../../../data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
freq_offset_file = '../output/Rev07E_RSS_2005_123_X43_E/freq_offset.txt'
geo_file = '../output/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB'
cal_file = '../output/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_CAL.TAB'
dlp_file = '../output/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_DLP.TAB'

geo_file = 'TEST_GEO.TAB'
cal_file = 'TEST_CAL.TAB'
dlp_file = 'TEST_DLP.TAB'

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

# GUI PARAMETERS USED:
#spm_include = 30250, 32600 ; 33520, 33890 ; 33990, 40200
#k = 5
#freespace_spm = 30130, 30800 ; 31782, 31790 ; 34074, 34083 ; 34240, 34252 ; 35275, 35295 ; 35530, 40000
#knots_spm = 30500, 31786, 34079, 34246, 35285, 37000
#k = 2

spm_include = [[30250, 32600], [33520, 33890], [33990, 40200]]
k_f_resid = 5
freespace_spm = [[30130, 30800], [31782, 31790], [34074, 34083],
    [34240, 34252], [35275, 35295], [35530, 40000]]
knots_spm = [30500, 31786, 34079, 34246, 35285, 37000]
k_power = 2
USE_GUI = True

# **END OF USER INPUT**

rsr_inst = RSRReader(rsr_file)

geo_inst = Geometry(rsr_inst, kernels)

f_spm, f_offset, freq_offset_history = calc_freq_offset(rsr_inst,
    freq_offset_file=freq_offset_file)
#f_offset_load = np.loadtxt(freq_offset_file)
#f_spm = f_offset_load[:, 0]
#f_offset = f_offset_load[:, 1]

fit_inst = FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset,
    f_uso, kernels, k=k_f_resid, spm_include=spm_include, USE_GUI=USE_GUI)
spm_vals, IQ_c = fit_inst.get_IQ_c()

norm_inst = Normalization(spm_vals, IQ_c, geo_inst, rsr_inst)
spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
    freespace_spm=freespace_spm, knots_spm=knots_spm, k=k_power,
    USE_GUI=USE_GUI)

cal_inst = Calibration(fit_inst, norm_inst, geo_inst)

make_geo_file(geo_inst, geo_file)
make_cal_file(cal_inst, cal_file)

# Beginning of quick-look steps - from here on out, you can run these steps
# without running the above, since the geo and cal files take care of the
# above steps

#cal_inst = MakeCalInst(cal_file)
norm_diff_inst = NormDiff(rsr_inst, dr_km_desired, geo_inst, cal_inst)
make_dlp_file(norm_diff_inst, dlp_file)
pdb.set_trace()
