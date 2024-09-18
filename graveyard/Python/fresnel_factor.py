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
#       Computes the Fresnel factor that occurs in MTR86.                      #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
import numpy

def fresnel_factor(rho, rho0, opening, azimuth):
    """
        Function: psi_factor_fast
        Purpose:  Calculate the first iteration of Newton-Raphson for psi with
            respect to phi using previously calculated sines and cosines.
        Variables:
            r:    Ring radius variable.
            r0:   Ring intercept point.
            b:    Rring opening angle.
            phi0: Ring azimuth angle.
    """
    cosb = numpy.cos(opening)
    sinphi = numpy.sin(azimuth)
    cosphi = numpy.cos(azimuth)
    cosb_sq = cosb*cosb
    sinphi_sq = sinphi*sinphi
    numerator = cosb_sq * cosphi * sinphi * (rho - rho0)
    denominator = (1.0 - cosb_sq * sinphi_sq) * rho0
    return numerator / denominator
