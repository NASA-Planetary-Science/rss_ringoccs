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
#       Computes the diffraction pattern through an inverted square well.      #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
from .square_well_diffraction import square_well_diffraction

def inverted_square_well_diffraction(x_val, left, right, fresnel_scale):
    """
        Function:
            inverted_square_well_diffraction
        Purpose:
            Computes diffraction pattern through an inverted square well.
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
    return 1.0 - square_well_diffraction(x_val, left, right, fresnel_scale)
