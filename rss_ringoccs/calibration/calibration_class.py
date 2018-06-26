"""

calibration_class.py

Purpose: Class for calibration parameters. This doesn't really do anything. It
         just puts everything calibration-related into one place

Revisions:
        calibration_class.py
    2018 Jun 11 - gsteranka - Original version
"""


import numpy as np
import os
import platform
from scipy.interpolate import interp1d
import spiceypy as spice
import time


class Calibration(object):
    """
    Purpose:
    Clump together everything calibration-related into a calibration instance.
    Use this if you haven't made a calibration file yet. If you have made a
    calibration file, use the "MakeCalInst" class instead.

    Attributes:
        t_oet_spm_vals (np.ndarray): SPM values from calibration file
        f_sky_hz_vals (np.ndarray): Predicted sky frequency values from
            calibration file
        f_sky_resid_fit_vals (np.ndarray): Fit to residual frequency from
            calibration file
        p_free_vals (np.ndarray): Freespace power spline fit values from
            calibration file
        f_offset_fit_vals (np.ndarray): Fit to frequency offset from
            calibration files
        history (dict): Dictionary with information of the run
    """

    def __init__(self, fit_inst, norm_inst, geo_inst, dt_cal=1.0):
        """
        Args:
            fit_inst: Instance of the FreqOffsetFit class
            norm_inst: Instance of the Normalization class
            geo_inst: Instance of the Geometry class
            dt_cal (float): Desired time Spacing between points
        """

        spm_geo = np.asarray(geo_inst.t_oet_spm_vals)
        rho_km_geo = np.asarray(geo_inst.rho_km_vals)
        phi_rad_vals = np.asarray(geo_inst.phi_ora_deg_vals)*spice.rpd()
        B_rad_vals = np.asarray(geo_inst.B_deg_vals)*spice.rpd()
        D_km_vals = np.asarray(geo_inst.D_km_vals)
        rho_dot_kms_vals = np.asarray(geo_inst.rho_dot_kms_vals)
        F_km_vals = np.asarray(geo_inst.F_km_vals)

        f_spm, f_sky_pred = fit_inst.get_f_sky_pred()
        f_spm, f_sky_resid_fit = fit_inst.get_f_sky_resid_fit()
        f_spm, f_offset_fit = fit_inst.get_f_offset_fit()

        # SPM for calibration parameters
        n_pts_cal = round((spm_geo[-1] - spm_geo[0])/dt_cal) + 1
        spm_cal = spm_geo[0] + dt_cal*np.arange(n_pts_cal)

        # Evaluate f_sky_pred at spm_cal
        f_sky_pred_func = interp1d(f_spm, f_sky_pred, fill_value='extrapolate')
        f_sky_pred_cal = f_sky_pred_func(spm_cal)

        # Evaluate f_sky_resid_fit at spm_cal
        f_sky_resid_fit_func = interp1d(f_spm, f_sky_resid_fit,
            fill_value='extrapolate')
        f_sky_resid_fit_cal = f_sky_resid_fit_func(spm_cal)

        # Evaluate f_offset_fit at spm_cal
        f_offset_fit_func = interp1d(f_spm, f_offset_fit,
            fill_value='extrapolate')
        f_offset_fit_cal = f_offset_fit_func(spm_cal)

        # Evaluate spline fit at spm_cal. Assumes you already made a
        #     satisfactory spline fit
        dummy_spm, p_free_cal = norm_inst.get_spline_fit(spm_fit=spm_cal,
        USE_GUI=False)

        self.t_oet_spm_vals = spm_cal
        self.f_sky_hz_vals = f_sky_pred_cal
        self.f_sky_resid_fit_vals = f_sky_resid_fit_cal
        self.p_free_vals = p_free_cal
        self.f_offset_fit_vals = f_offset_fit_cal
        self.__set_history(fit_inst, norm_inst, geo_inst, dt_cal)


    def __set_history(self, fit_inst, norm_inst, geo_inst, dt_cal):

        input_var_dict = {'fit_inst': fit_inst.history,
            'norm_inst': norm_inst.history, 'geo_inst': geo_inst.history}
        input_kw_dict = {'dt_cal': dt_cal}

        hist_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__.split('/')[-1],
            'Source Directory': __file__.rsplit('/',1)[0] +'/',
            'Input Variables': input_var_dict,
            'Input Keywords': input_kw_dict}

        self.history = hist_dict
