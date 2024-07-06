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
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# Ignore this very silly warning.
# pylint: disable = too-many-arguments
import numpy
from .resolution_inverse import resolution_inverse

def window_width_and_prange(res, normeq, fsky, fres, rho_dot, sigma, bfac=True):
    """
        Purpose:
            Compute the window width as a function of ring radius.
            This is given from MTR86 Equations 19, 32, and 33.
        Variables:
            :res (*float*):
                The requested resolution.
            :normeq (*float*):
                The normalized equivalent width. Unitless.
            :fsky (*float* or *numpy.ndarray*):
                The sky frequency.
            :fres (*float* or *numpy.ndarray*):
                The Fresnel scale.
            :rdot (*float*) or (*numpy.ndarray*):
                The time derivative of the ring radius.
        Output:
            :w_vals (*numpy.ndarray*):
                The window width as a function of ring radius.
    """
    if bfac:
        omega = 2.0 * numpy.pi * fsky
        alpha = numpy.square(omega * sigma) / (2.0 * rho_dot)
        p_vals = res / (alpha * numpy.square(fres))

        # Create a variable specifying where P>1 occurs.
        p_range = (p_vals > 1.0).nonzero()[0]

        if numpy.size(p_range) == 0:
            raise IndexError(
                """
                \r\tError Encountered: rss_ringoccs
                \r\t\tdiffrec.special_functions.window_width\n
                \r\tThe P parameter in window width computation is less than
                \r\tone for the entirety of the data set. Either
                \r\trho_dot_km_vals is too small, tor F_km_vals is too large.
                \r\tRequest a coarser resolution, or check your data for
                \r\terrors. Alternatively, you may set bfac=False as a keyword.
                \r\tThis may result in inaccurate window width.
                """
            )

        p_vals = p_vals[p_range]
        alpha = alpha[p_range]

        w_vals = numpy.zeros(numpy.size(rho_dot))
        w_vals[p_range] = resolution_inverse(p_vals) / alpha

    else:
        w_vals = 2.0*numpy.square(fres)/res
        p_range = (fres > 0.0).nonzero()[0]

    w_vals *= normeq

    return w_vals, p_range
