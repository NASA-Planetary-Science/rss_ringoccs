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
#       Computes the diffraction pattern through a square well.                #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
from .fresnel_cos import fresnel_cos
from .fresnel_sin import fresnel_sin

def square_well_diffraction(x_val, left, right, fresnel_scale):
    """
        Function:
            square_well_diffraction
        Purpose:
            Computes the solution of diffraction through a square well.
        Arguments:
            x_val:
                The independent variable.
            left:
                The left-most endpoint of the square well.
            right:
                The right-most endpoint of the square well.
            fresnel_scale:
                The Fresnel scale.
        Output:
            square_well:
                Diffraction pattern of a square well.
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    factor = 0.5 - 0.5j
    rcpr_f_val = 1.0 / fresnel_scale
    left_factor = (left - x_val) * rcpr_f_val
    right_factor = (right - x_val) * rcpr_f_val

    fc_factor = fresnel_cos(right_factor) - fresnel_cos(left_factor)
    fs_factor = fresnel_sin(right_factor) - fresnel_sin(left_factor)

    return factor * (fc_factor + 1.0j*fs_factor)
