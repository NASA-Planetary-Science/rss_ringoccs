"""

make_cal_inst.py

Purpose: Create an instance for calibration parameters. Just for NormDiff class,
         for the purpose of making all inputs an instance of something. This is
         a replacement for the "Calibration" class in calibration_class.py if
         you have already made a calibration file

Revisions:
        make_cal_inst.py
    2018 Jun 11 - gsteranka - Original version
"""

import numpy as np
import os
import platform
import time


class MakeCalInst(object):
    """
    Purpose:
    Make a calibration instance using a calibration file

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

    Notes:
        [1] Works on 5-COLUMN CAL FILES ONLY!!! Need the 5th column of
            frequency offset fit to avoid running whole fit routine again
    """

    def __init__(self, cal_file):
        """
        Args:
            cal_file (str): Full path name of calibration file
        """

        cal = np.loadtxt(cal_file, delimiter=',')

        t_oet_spm_vals = cal[:, 0]
        f_sky_pred_vals = cal[:, 1]
        f_sky_resid_fit_vals = cal[:, 2]
        p_free_vals = cal[:, 3]
        f_offset_fit_vals = cal[:, 4]

        self.t_oet_spm_vals = t_oet_spm_vals
        self.f_sky_hz_vals = f_sky_pred_vals
        self.f_sky_resid_fit_vals = f_sky_resid_fit_vals
        self.p_free_vals = p_free_vals
        self.f_offset_fit_vals = f_offset_fit_vals
        self.__set_history(cal_file)


    def __set_history(self, cal_file):

        input_var_dict = {'cal_file': cal_file}

        hist_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__,
            'Input Variables': input_var_dict}
        self.history = hist_dict
