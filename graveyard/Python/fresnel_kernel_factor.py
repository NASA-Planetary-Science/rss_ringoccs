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
#       Computes the Fresnel kernel factor from MTR86.                         #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy

def fresnel_kernel_factor(rho, rho0, opening, phi0):
    """
        Function:
            fresnel_kernel_factor
        Purpose: 
            Calculate the first iteration of Newton-Raphson for
            psi with respect to phi using previously calculated
            sines and cosines.
        Variables:
            rho:
                Ring radius variable.
            rho0:
                Ring intercept point.
            opening:
                Ring opening angle.
            phi0:
                Ring azimuth angle.
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2024/07/06: Ryan Maguire
                Cleaned up and added to graveyard.
    """
    cb_sq = numpy.square(numpy.cos(opening))
    cp0 = numpy.cos(phi0)
    sp0 = numpy.sin(phi0)
    sp0_sq = numpy.square(sp0)

    return (cb_sq*cp0*sp0 / (1.0 - cb_sq * sp0_sq)) * (rho - rho0) / rho0
