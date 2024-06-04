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
#       Gets the pre-computed normalized equivalent width from a string.       #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
from .windows import FUNCTIONS

def get_norm_eq(wtype):
    """
        Function:
            get_norm_eq
        Purpose:
            Compute the Normalized Equivalent Width from a given WTYPE.
        Variables:
            wtype:
                The name of the window used for processing.
        Output:
            norm_eq:
                The normalized equivalent width of wtype.
        Notes:
            [1] The normalized equivalent width is pre-computed using
                a window with infinite resolution. That is, the
                intergral is perform on the functions exactly, not
                approximated with Riemann sums. This value is more
                accurate than the method used in the compute_norm_eq
                function, however only a few select windows are
                included. These windows are:
                    [1] Rectangular Window..................'rect'
                    [2] Squared Cosine Window...............'coss'
                    [3] Kaiser-Bessel 2.5 Window............'kb25'
                    [4] Kaiser-Bessel 3.5 Window............'kb35'
                    [5] Modified Kaiser-Bessel 2.5 Window...'kbmd'
            [2] The numerical values that this function returns are
                slightly different than those quoted in the reference
                below. The MTR86 paper only evaluates to two decimals
                whereas we have done double precision to 8 decimals.
            [3] The input must be a string. The function is neither
                case nor space sensitive.
            [4] The normalized equivalent width is a unitless value.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
    """
    if not isinstance(wtype, str):
        raise TypeError("Input must be a string.")

    # Get rid of any whitespace and quotations marks.
    window_type = wtype.replace(" ", "").replace("'", "").replace("\"", "")

    # Make the string lower case.
    window_type = window_type.lower()

    if window_type not in FUNCTIONS:
        raise ValueError("Invalid window type.")

    return FUNCTIONS[window_type]["normeq"]
