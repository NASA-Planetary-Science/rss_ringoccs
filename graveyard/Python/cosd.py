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
#       Computes cosine in degrees.                                            #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy

def cosd(x_val):
    """
        Function:
            cosd
        Purpose:
            Computes cosine in degrees.
        Arguments:
            x_val:
                The independent variable, in degrees.
        Output:
            cos_x:
                The cosine of the input.
    """
    # Converts degrees to radians.
    conversion_factor = 0.017453292519943295

    # Reduce argument to lie between -360 and 360.
    # Multiplying the conversion factor by a very large
    # number can lead to precision loss. This avoids
    # such problems.
    reduced = numpy.fmod(x_val, 360.0)

    # Convert to radians and compute.
    return numpy.cos(reduced * conversion_factor)
