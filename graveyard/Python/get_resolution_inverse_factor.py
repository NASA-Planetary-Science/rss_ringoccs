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
#       Computes the inverse factor for the resolution function.               #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy
from .resolution_inverse import resolution_inverse
from .window_range_data import WindowRangeData

def get_resolution_inverse_factor(window_data):
    """
        Function:
            get_resolution_inverse_factor
        Purpose:
            Computes the inverse factor for the window function.
        Arguments:
            window_data:
                Instance of the WindowRangeData class.
        Output:
            inverse_factor:
                The factor used in the inverse function for
                the resolution as a function of window width.
    """
    if not isinstance(window_data, WindowRangeData):
        raise TypeError("Input not an instance of the WindowRangeData class.")

    omega = 2.0 * numpy.pi * window_data.fsky
    denom = 2.0 * window_data.rho_dot
    numer = numpy.square(omega * window_data.sigma)
    rcpr_alpha = denom / numer
    fres_sq = numpy.square(window_data.fres)

    factor = window_data.res * rcpr_alpha / fres_sq

    # The inverse exists only if factor > 1.
    if numpy.min(factor) < 1.0001:
        raise ValueError(
            "Bad Points found in requested range.\n"
            "Either rho_dot_kms_vals, F_km_vals, or res_km is to small.\n"
            "Exclude points or request coarser resolution.\n"
        )

    return numpy.abs(resolution_inverse(factor) * rcpr_alpha)
