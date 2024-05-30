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
#       Extracts reconstructed data from a CSV file.                           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# Pylint will complain about "duplicate code" because the DLP and TAU files
# have 5 lines in common, in the same order, for their "NAMES" variable.
# This cannot be helped, this is how the data is structured.
# Silence this warning.
# pylint: disable = duplicate-code
import pandas

# The names for the columns in the order that they appear.
TAU_NAMES = [
    "rho_km_vals",
    "rho_km_pole_corr_vals",
    "rho_km_offsett_vals",
    "phi_rl_deg_vals",
    "phi_ora_deg_vals",
    "raw_tau_vals",
    "phase_deg_vals",
    "raw_tau_threshold_vals",
    "spm_vals",
    "t_ret_spm_vals",
    "t_set_spm_vals",
    "B_deg_vals"
]

def get_tau(tau, verbose = True):
    """
        Function:
            get_tau
        Purpose:
            Extracts reconstructed data from a CSV file.
        Arguments:
            tau (str):
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
    if not isinstance(tau, str):
        raise TypeError("tau must be a string: '/path/to/tau'")

    if not isinstance(verbose, bool):
        raise TypeError("verbose must be Boolean: True/False")

    if verbose:
        print("\tExtracting Tau Data...")

    csv_data = pandas.read_csv(tau, delimiter = ',', names = TAU_NAMES)

    if verbose:
        print("\tTau Data Complete.")

    return csv_data
