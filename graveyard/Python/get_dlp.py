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
#       Extracts DLP data from a CSV file.                                     #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import pandas

# The names for the columns in the order that they appear.
DLP_NAMES = [
    "rho_km_vals",
    "rho_corr_pole_km_vals",
    "rho_corr_timing_km_vals",
    "phi_rl_deg_vals",
    "phi_ora_deg_vals",
    "raw_tau_vals",
    "phase_deg_vals",
    "raw_tau_threshold_vals",
    "t_oet_spm_vals",
    "t_ret_spm_vals",
    "t_set_spm_vals",
    "B_deg_vals"
]

def get_dlp(dlp, verbose = True):
    """
        Function:
            get_dlp
        Purpose:
            Extracts DLP data from a CSV file.
        Arguments:
            dlp (str):
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
    if not isinstance(dlp, str):
        raise TypeError("dlp must be a string: '/path/to/dlp'")

    if not isinstance(verbose, bool):
        raise TypeError("verbose must be Boolean: True/False")

    if verbose:
        print("\tExtracting DLP Data...")

    csv_data = pandas.read_csv(dlp, delimiter = ',', names = DLP_NAMES)

    if verbose:
        print("\tDLP Data Complete.")

    return csv_data
