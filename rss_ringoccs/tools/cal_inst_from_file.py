"""

cal_inst_from_file.py

Purpose: Create an instance for calibration parameters. Just for NormDiff class,
         for the purpose of making all inputs an instance of something. This is
         a replacement for the "Calibration" class in calibration_class.py if
         you have already made a calibration file

WHERE TO GET NECESSARY INPUT:
    cal_file: Output from the PDS3 CAL file writer. This isn't currently on
        GitHub, but you can find it under the TC2017 directory at
        TC2017/jfong/programs/occgeo/jwf_pds3_cal_series_v2.py
    rsr_inst: Use an instance of the RSRReader class, found inside of
        rss_ringoccs/rsr_reader/rsr_reader.py

Revisions:
        make_cal_inst.py
    2018 Jun 11 - gsteranka - Original version
    2018 Jun 27 - gsteranka - Adjust to sky frequency in CAL has frequency
                              offset fit added on top of it
    2018 Aug 10 - jfong - rename file to cal_inst_from_file.py from
                          make_cal_inst.py
                        - rename class to CreateCalInst from MakeCalInst
"""

import numpy as np
import os
import platform
import sys
import time

sys.path.append('../..')
from rss_ringoccs.rsr_reader.rsr_reader import RSRReader
sys.path.remove('../..')


class CreateCalInst(object):
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
    """

    def __init__(self, cal_file, rsr_inst):
        """
        Args:
            cal_file (str): Full path name of calibration file
            rsr_inst: Instace of the RSRReader class
        """

        try:
            cal = np.loadtxt(cal_file, delimiter=',')
        except FileNotFoundError:
            sys.exit('ERROR (MakeCalInst): File not found')

        if type(rsr_inst) != RSRReader:
            sys.exit('ERROR (MakeCalInst): rsr_inst input needs to be an '
                + 'instance of the RSRReader class')

        t_oet_spm_vals = cal[:, 0]
        f_sky_hz_vals = cal[:, 1]
        f_sky_resid_fit_vals = cal[:, 2]
        p_free_vals = cal[:, 3]

        dummy_spm, f_sky_pred_file_vals = rsr_inst.get_f_sky_pred(
            f_spm=t_oet_spm_vals)
        f_offset_fit_vals = f_sky_hz_vals - f_sky_pred_file_vals

        self.t_oet_spm_vals = t_oet_spm_vals
        self.f_sky_hz_vals = f_sky_hz_vals
        self.f_sky_resid_fit_vals = f_sky_resid_fit_vals
        self.p_free_vals = p_free_vals
        self.f_offset_fit_vals = f_offset_fit_vals
        self.__set_history(cal_file)


    def __set_history(self, cal_file):
        """
        Set history attribute
        """

        input_var_dict = {'cal_file': cal_file}

        hist_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__,
            'Input Variables': input_var_dict}
        self.history = hist_dict
