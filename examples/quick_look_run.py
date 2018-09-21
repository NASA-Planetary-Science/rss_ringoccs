"""

quick_look_run.py

Purpose: Example rss_ringoccs/ 'quick-look' process script.

Revisions:
    2018 Sep 19 - jfong - original
"""

import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import pdb

# ***** Begin user input *****
data_dir = '/Volumes/jfong001/Research/TC2017/data/Cassini_RSS_Ring_Profiles_2018_Archive/Rev7E/Rev007E_RSS_2005_123_X43_E/'
geo_file = data_dir + 'RSS_2005_123_X43_E_GEO.TAB'
cal_file = data_dir + 'RSS_2005_123_X43_E_CAL.TAB'
dlp_file = data_dir + 'RSS_2005_123_X43_E_DLP_500M.TAB'

verbose = True
res_km = 1.0
inversion_range = 'all'

# ***** End user input *****

dlp_inst = rss.tools.ExtractCSVData(geo_file, cal_file, dlp_file,
        verbose=verbose)

tau_inst = rss.diffcorr.DiffractionCorrection(dlp_inst, res_km,
        rng=inversion_range, verbose=verbose)

pdb.set_trace()

