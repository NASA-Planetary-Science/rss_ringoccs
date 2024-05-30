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
#       Extracts calibration data from a CSV file.                             #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import pandas

# The names for the columns in the order that they appear.
CAL_NAMES = [
    "spm_vals",
    "f_sky_pred_vals",
    "f_sky_resid_fit_vals",
    "p_free_vals"
]

def get_cal(cal, verbose = True):
    """
        Function:
            get_cal
        Purpose:
            Extracts calibration data from a CSV file.
        Arguments:
            cal (str):
                The path to the CSV file.
        Keywords:
            verbose (bool):
                Boolean for printing out messages.
        Output:
            csv_data (pandas.core.frame.DataFrame):
                The data in the CSV file.
        History:
            2018: Ryan Maguire
                Original draft.
            2024/05/30: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    if not isinstance(cal, str):
        raise TypeError("cal must be a string: '/path/to/cal'")

    if not isinstance(verbose, bool):
        raise TypeError("verbose must be Boolean: True/False")

    if verbose:
        print("\tExtracting CAL Data...")

    csv_data = pandas.read_csv(cal, delimiter = ',', names = CAL_NAMES)

    if verbose:
        print("\tCAL Data Complete.")

    return csv_data
