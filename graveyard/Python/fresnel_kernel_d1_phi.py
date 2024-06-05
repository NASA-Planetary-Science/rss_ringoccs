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
#       Computes d Psi / d Phi, the derivative of the Fresnel kernel.          #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
# Ignore silly Pylint warnings.
# pylint: disable = too-many-arguments
# pylint: disable = too-many-locals
import numpy

def fresnel_kernel_d1_phi(rho, rho0, phi, phi0, opening, dist):
    """
        Function:
            fresnel_kernel_d1_phi
        Purpose:
            Calculate the partial derivative of the Fresnel
            kernel with respect to the azimuth angle.
        Variables:
            rho:
                Ring radius variable.
            rho0:
                Ring intercept point.
            phi:
                The ring azimuth angle for rho.
            phi0
                The ring azimuth angle for rho0.
            opening:
                The ring opening angle.
            dist:
                RIP-Spacecraft distance.
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2024/06/05: Ryan Maguire
                Cleaned up and added to graveyard.
    """
    # Precompute trig values.
    cosb = numpy.cos(opening)
    sinp = numpy.sin(phi)
    cosp = numpy.cos(phi)
    sinp0 = numpy.sin(phi0)
    cosp0 = numpy.cos(phi0)

    # Precompute divisions and products.
    rcpr_d = 1.0 / dist
    rcpr_d_sq = numpy.square(rcpr_d)
    rho_prod = 2.0*rho*rho0
    eta_term = -rho_prod * (sinp*sinp0 + cosp*cosp0)

    # Compute the derivative.
    xi_factor = cosb * rcpr_d * (rho0*cosp0 - rho*cosp)
    eta_factor = (numpy.square(rho0) + numpy.square(rho) + eta_term)*rcpr_d_sq
    deriv1 = rho*cosb*sinp*rcpr_d
    deriv2 = rho_prod*(sinp*cosp0 - cosp*sinp0)*rcpr_d_sq
    numer = deriv1 + 0.5*deriv2
    denom = numpy.sqrt(1.0 + 2.0*xi_factor + eta_factor)
    return -deriv1 + numer / denom
