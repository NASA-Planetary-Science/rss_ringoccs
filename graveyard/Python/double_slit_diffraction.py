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
#       Creates the double slit diffraction pattern.                           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy

def double_slit_diffraction(x_val, dist, width, separation):
    """
        Function:
            double_slit_diffraction
        Purpose:
            Function for the Fraunhofer double slit diffraction formula.
        Arguments:
            x_val (float):
                Location on the x-axis.
            dist (float):
                Distance on the y-axis from the observer to the slit.
            width (float):
                The width of the slit.
            separation (float):
                Distance between slits.
        Output:
            diffraction (float):
                The diffraction pattern using Fraunhofer optics.
        History:
            2018:   Ryan Maguire
                Original draft.
            2024/05/30: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    first_arg = x_val * width / dist
    second_arg = numpy.pi * separation * x_val / dist

    first_term = numpy.square(numpy.sinc(first_arg))
    second_term = numpy.square(numpy.sin(2.0 * second_arg))
    third_term = 4.0 * numpy.square(numpy.sin(second_term))

    return first_term * second_term / third_term
