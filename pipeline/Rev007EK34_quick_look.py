"""

Rev007EK34_quick_look.py

Purpose: Demonstrate the quick-look procedure on a set of pre-made data files
         made from the Rev007EK34_e2e.py script

Revisions:
        Rev007EK34_quick_look.py
    2018 Aug 28 - gsteranka - Original version
"""

import matplotlib.pyplot as plt
import sys

sys.path.append('..')
import rss_ringoccs as rss
sys.path.remove('..')

# Replace geo_file, cal_file, and dlp_file with the name of your output
#     from Rev007EK34_e2e.py
output_directory = '../output/rev7E_K34_e2e_output/'
geo_file = output_directory + 'RSS_2005_123_K34_E_GEO_20180828.TAB'
cal_file = output_directory + 'RSS_2005_123_K34_E_CAL_20180828.TAB'
dlp_file = output_directory + 'RSS_2005_123_K34_E_DLP_20180828.TAB'

res_km = 1.0
inversion_range = [117700, 117900]

csv_data = rss.tools.ExtractCSVData(geo_file, cal_file, dlp_file)

tau_inst = rss.diffcorr.DiffractionCorrection(csv_data, res_km,
    rng=inversion_range)

plt.plot(tau_inst.rho_km_vals, tau_inst.power_vals)
plt.title('Rev007E K34 Huygens Ringlet')
plt.xlabel('Rho (km)')
plt.ylabel('Power')
plt.show()
