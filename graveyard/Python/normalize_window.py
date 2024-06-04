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
#       Normalizes a window by the free-space integral.                        #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy
from .check_pos_real import check_pos_real

def normalize_window(bin_width, ker, f_scale, check_error = True):
    """
        Function:
            normalize
        Purpose:
            Normalizes a window function by the window width.
    """

    if check_error:
        if not check_pos_real(bin_width):
            raise ValueError("bin_width must be a positive real number.")

        if not check_pos_real(f_scale):
            raise ValueError("f_scale must be a positive real number")

        if numpy.size(ker) == 0:
            raise TypeError("ker must be a non-empty array.")

    # Freespace Integral.
    denom = numpy.abs(numpy.sum(ker) * bin_width)

    # Normalization factor from the Fresnel scale.
    numer = 1.4142135623730951 * f_scale

    # Normalization factor.
    return numer / denom
