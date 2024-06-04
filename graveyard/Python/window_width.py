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
#       Computes the required window width as a function of the resolution.    #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy
from .get_resolution_inverse_factor import get_resolution_inverse_factor
from .window_range_data import WindowRangeData

def window_width(window_data):
    """
        Function:
            window_width
        Purpose:
            Compute the window width as a function of
            resolution and ring radius in the plane.
        Arguments:
            window_data:
                An instance of the WindowRangeData class.
        Output:
            w_vals:
                The window width as a function of ring radius.
    """

    if not isinstance(window_data, WindowRangeData):
        raise TypeError("Input not an instance of the WindowRangeData class.")

    # Check if the b-factor is to be used.
    if window_data.bfac:
        inverse_factor = get_resolution_inverse_factor(window_data)
        w_vals = window_data.normeq * inverse_factor

    # Without the b-factor, the formula is very simple.
    else:
        fres_sq = numpy.square(window_data.fres)
        w_vals = 2.0 * window_data.normeq * fres_sq / window_data.res

    return w_vals
