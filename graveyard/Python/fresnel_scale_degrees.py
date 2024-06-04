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
#       Computes the Fresnel scale.                                            #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/04/14                                                         #
################################################################################
"""
import numpy
from .cosd import cosd
from .sind import sind

def fresnel_scale_degrees(wavelength, distance, azimuth_angle, opening_angle):
    """
        Function:
            fresnel_scale_degrees
        Purpose:
            Compute the Fresnel Scale from lambda, D, Phi, and B, in degrees.
        Variables:
            wavelength:
                Wavelength of the incoming signal.
            distance:
                ring-to-spacecraft distance.
            azimuth_angle:
                Ring azimuth angle.
            opening_angle:
                Ring opening angle.
        Output:
            fresnel_scale:
                The Fresnel scale.
        History:
            2018/04/14: Ryan Maguire
                Translated from IDL.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    cos_b = cosd(opening_angle)
    sin_b = sind(opening_angle)
    sin_p = sind(azimuth_angle)
    numer = 1.0 - numpy.square(cos_b * sin_p)
    denom = numpy.square(sin_b)
    return numpy.sqrt(0.5 * wavelength * distance * numer / denom)
