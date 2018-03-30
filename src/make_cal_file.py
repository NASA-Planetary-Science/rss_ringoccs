"""

make_cal_and_obs_file.py

Purpose: Write out calibration ("cal") and observations ("obs") files for
         the occultation. These will be used in the "Quick-Look" process,
         and produced in the "End-to-End" process

NOTE: Dependent on final format of geometry instances ("geo_inst")

Revisions:
      gjs_make_cal_and_obs_file.py
   Mar 19 2018 - gsteranka - Original version
      make_cal_and_obs_file.py
   Mar 20 2018 - gsteranka - Copied to current version and deleted debug steps
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import spiceypy as spice

def make_cal_file(cal_file, fit_inst, norm_inst, geo_inst, dt_cal=1.0):
    """Save a cal file from extracted data.

    Args:
        cal_file (str): Full path name of cal file to be made
        fit_inst: Instance of the FreqOffsetFit class in freq_offset_fit.py
        norm_inst: Instance of the Normalization class in power_normalization.py
        geo_inst: Mock instance created from a geometry file using
            geo_file_into_instance.py
        dt_cal (float): Time spacing of the calibration file. Default is 1 sec
    """

    # TODO (gsteranka): Update attribute names once I get geometry routines
    spm_geo = np.asarray(geo_inst.t_oet_spm_vals)
    rho_km_geo = np.asarray(geo_inst.rho_km_vals)
    phi_rad_vals = np.asarray(geo_inst.phi_ora_deg_vals)*spice.rpd()
    B_rad_vals = np.asarray(geo_inst.B_deg_vals)*spice.rpd()
    D_km_vals = np.asarray(geo_inst.D_km_vals)
    rho_dot_kms_vals = np.asarray(geo_inst.rho_dot_kms_vals)
    F_km_vals = np.asarray(geo_inst.F_km_vals)

    # For cal file
    (f_spm, f_sky_pred) = fit_inst.get_f_sky_pred()
    (f_spm, f_sky_resid_fit) = fit_inst.get_f_sky_resid_fit()
    (f_spm, f_offset_fit) = fit_inst.get_f_offset_fit()

    # SPM for cal file
    n_pts_cal = round((spm_geo[-1] - spm_geo[0])/dt_cal) + 1
    spm_cal = spm_geo[0] + dt_cal*np.arange(n_pts_cal)

    # Evaluate f_sky_pred at SPM_cal
    f_sky_pred_func = interp1d(f_spm, f_sky_pred, fill_value='extrapolate')
    f_sky_pred_cal = f_sky_pred_func(spm_cal)

    # Evaluate f_sky_resid_fit at SPM_cal
    f_sky_resid_fit_func = interp1d(f_spm, f_sky_resid_fit,
                                    fill_value='extrapolate')
    f_sky_resid_fit_cal = f_sky_resid_fit_func(spm_cal)

    # Evaluate f_offset_fit at spm_cal
    f_offset_fit_func = interp1d(f_spm, f_offset_fit, fill_value='extrapolate')
    f_offset_fit_cal = f_offset_fit_func(spm_cal)

    # Evaluate spline fit at spm_cal
    p_free_cal = norm_inst.get_spline_fit(spm_cal)

    np.savetxt(cal_file, np.c_[spm_cal, f_sky_pred_cal, f_sky_resid_fit_cal,
                               p_free_cal, f_offset_fit_cal],
               fmt='%32.16f '*5)

    return None

if __name__ == '__main__':
    pass
