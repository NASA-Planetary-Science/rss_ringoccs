"""

make_cal_file.py

Purpose: Write out calibration ("cal") and observations ("obs") files for
         the occultation. These will be used in the "Quick-Look" process,
         and produced in the "End-to-End" process

NOTE: Dependent on final format of geometry instances ("geo_inst")

Revisions:
      gjs_make_cal_and_obs_file.py
   2018 Mar 19 - gsteranka - Original version
      make_cal_file.py
   2018 Mar 20 - gsteranka - Copied to current version and deleted debug steps
   2018 May 23 - gsteranka - Edited to save CSV file, rather than white-space
                             delimited
   2018 Jun 11 - gsteranka - Accept instance of the Calibration class in
                             calibration_class.py
"""

import numpy as np
from scipy.interpolate import interp1d
import spiceypy as spice

def make_cal_file(cal_inst, cal_file):
    """
    Purpose:
    Save a cal file from extracted data

    Args:
        cal_inst: Instance of Calibration class
        cal_file (str): Full path name of cal file to be made
    """

    spm_cal = cal_inst.t_oet_spm_vals
    f_sky_pred_cal = cal_inst.f_sky_pred_vals
    f_sky_resid_fit_cal = cal_inst.f_sky_resid_fit_vals
    p_free_cal = cal_inst.p_free_vals
    f_offset_fit_cal = cal_inst.f_offset_fit_vals

    np.savetxt(cal_file, np.c_[spm_cal, f_sky_pred_cal, f_sky_resid_fit_cal,
        p_free_cal, f_offset_fit_cal], fmt='%32.16f, '*4 + '%32.16f')

    return None

if __name__ == '__main__':
    pass
