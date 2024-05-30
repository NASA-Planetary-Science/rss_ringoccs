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
#       Creates the single slit diffraction pattern.                           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy

def single_slit_diffraction(x_val, dist, width):
    """
        Function:
            single_slit_diffraction
        Purpose:
            Function for the Fraunhofer single slit diffraction formula.
        Arguments:
            x_val (float):
                Location on the x-axis.
            dist (float):
                Distance on the y-axis from the observer to the slit.
            width (float):
                The width of the slit.
        Output:
            diffraction (float):
                The diffraction pattern using Fraunhofer optics.
        History:
            2018:   Ryan Maguire
                Original draft.
            2024/05/30: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    arg = x_val * width / dist
    return numpy.square(numpy.sinc(arg))
