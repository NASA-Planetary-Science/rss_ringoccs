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
#       Computes the Fresnel kernel.                                           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# Ignore silly Pylint warnings.
# pylint: disable = too-many-arguments
# pylint: disable = too-many-locals
import numpy

def fresnel_kernel(kd_val, rho, rho0, phi, phi0, opening, dist):
    """
        Function:
            fresnel_kernel
        Purpose:
            Calculate the Fresnel kernel.
        Variables:
            kd_val:
                The wavenumber scaled by the distance.
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
            2024/07/06: Ryan Maguire
                Cleaned up and added to graveyard.
    """
    # Precompute trig values.
    cosb = numpy.cos(opening)
    sinp = numpy.sin(phi)
    cosp = numpy.cos(phi)
    sinp0 = numpy.sin(phi0)
    cosp0 = numpy.cos(phi0)
    cos_p_p0 = cosp*cosp0 + sinp*sinp0

    # Precompute divisions and products.
    rcpr_d = 1.0 / dist
    rcpr_d_sq = numpy.square(rcpr_d)

    # Compute the Fresnel kernel.
    xi_factor = (rho * cosp - rho0 * cosp0) * cosb * rcpr_d
    eta_factor = (rho*rho + rho0*rho0 - 2.0*rho*rho0*cos_p_p0) * rcpr_d_sq
    psi0 = numpy.sqrt(1.0 - 2.0*xi_factor + eta_factor)
    return kd_val * (psi0 - 1 + xi_factor)
