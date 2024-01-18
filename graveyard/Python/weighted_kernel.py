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
#       Computes the Fresnel kernel weighted by a window function.             #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import numpy

def weighted_kernel(window_array, fresnel_psi):
    """
        Function:
            weighted_kernel
        Purpose:
            Computes the weighted Fresnel kernel function.
        Variables:
            window_array:
                The data for the window function.
            fresnel_psi:
                The geometric quantity from the Fresnel kernel.
        Output:
            kernel:
                The weighted Fresnel kernel.
    """

    return window_array * numpy.exp(1j * fresnel_psi)
