"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Purpose:                                                                   #
#       Create an instance for calibration parameters. Just for NormDiff       #
#       class, for the purpose of making all inputs an instance of something.  #
#       This is a replacement for the "Calibration" class in                   #
#       calibration_class.py if you have already made a calibration file       #
#       WHERE TO GET NECESSARY INPUT:                                          #
#           cal_file:                                                          #
#               Output from the PDS3 CAL file writer. This isn't currently on  #
#               GitHub, but you can find it under the TC2017 directory at      #
#               TC2017/jfong/programs/occgeo/jwf_pds3_cal_series_v2.py         #
#           rsr_inst:                                                          #
#               Use an instance of the RSRReader class, found inside of        #
#               rss_ringoccs/rsr_reader/rsr_reader.py                          #
################################################################################
#   Author:     Glenn Steranka                                                 #
#   Date:       2018/06/11                                                     #
################################################################################
"""
# Disable silly pylint warning.
# pylint: disable = too-few-public-methods
import os
import platform
import time
import numpy

class CreateCalInst:
    """
    Purpose:
        Make a calibration instance using a calibration file
    Attributes:
        t_oet_spm_vals (numpy.ndarray): SPM values from calibration file
        f_sky_hz_vals (numpy.ndarray): Predicted sky frequency values from
            calibration file
        f_sky_resid_fit_vals (numpy.ndarray): Fit to residual frequency from
            calibration file
        p_free_vals (numpy.ndarray): Freespace power spline fit values from
            calibration file
        f_offset_fit_vals (numpy.ndarray): Fit to frequency offset from
            calibration files
        history (dict): Dictionary with information of the run
    """

    def __init__(self, cal_file, rsr_inst):
        """
        Args:
            cal_file (str):
                Full path name of calibration file
            rsr_inst:
                Instace of the RSRReader class
        """
        cal = numpy.loadtxt(cal_file, delimiter = ',')

        t_oet_spm_vals = cal[:, 0]
        f_sky_hz_vals = cal[:, 1]
        f_sky_resid_fit_vals = cal[:, 2]
        p_free_vals = cal[:, 3]

        dummy_spm, f_sky_pred_file_vals = rsr_inst.get_f_sky_pred(
            f_spm = t_oet_spm_vals
        )

        f_offset_fit_vals = f_sky_hz_vals - f_sky_pred_file_vals

        self.t_oet_spm_vals = t_oet_spm_vals
        self.f_sky_hz_vals = f_sky_hz_vals
        self.f_sky_resid_fit_vals = f_sky_resid_fit_vals
        self.p_free_vals = p_free_vals
        self.f_offset_fit_vals = f_offset_fit_vals
        self.set_history(cal_file)

    def set_history(self, cal_file):
        """
        Set history attribute
        """

        input_var_dict = {'cal_file': cal_file}

        hist_dict = {
            'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__,
            'Input Variables': input_var_dict
        }

        self.history = hist_dict
