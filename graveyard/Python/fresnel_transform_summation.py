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
#       Computes the numerical Fresnel transform.                              #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/04/15                                                         #
################################################################################
"""
import numpy

def fresnel_transform_summation(t_vals, kernel, displacement, fresnel_scale):
    """
        Function:
            fresnel_transform_summation
        Purpose:
            Compute the discrete Fresnel transform by summation.
        Arguments:
            t_vals:
                The input function (impulse data).
            kernel:
                The Fresnel kernel.
            displacement:
                The size of a bin (r[1] - r[0]).
            fresnel_scale:
                The Frensel scale.
        Output:
            t_hat:
                The forward model for the diffraction pattern.
        History:
            2018/04/15: Ryan Maguire
                Translated from IDL.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    factor = (0.5 - 0.5j) * displacement / fresnel_scale
    return numpy.sum(kernel * t_vals) * factor
