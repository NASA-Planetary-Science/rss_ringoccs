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
#       Provides the standard rectangular window function.                     #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
import numpy
from .check_pos_real import check_pos_real
from .number_of_window_points import number_of_window_points

def rect(w_in, bin_width, check_error = True):
    """
        Function:
            rect
        Purpose:
            Create the rectangular window function.
        Arguments:
            w_in (float):
                Window width.
            bin_width (float):
                Width of one point (or one bin).
        Outputs:
            w_func (numpy.ndarray):
                The rectungular window function of width w_in
                and spacing bin_width between points.
        Dependencies:
            [1] numpy
        Notes:
            The rectangular window function is the unit function.
            That is, it is equal to one across it's entire domain. It
            is used in the Fresnel Inversion and Fresnel Transform
            process to act as a 'hat' function, zeroing out a
            function outside of the domain of rect's definition, and
            having no effect within that domain.
        References:
            [1] https://en.wikipedia.org/wiki/Window_function
            [2] https://en.wikipedia.org/wiki/Rectangular_function
            [3] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2018/05/16: Ryan Maguire
                Made variables lower case.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    if check_error:
        if (not check_pos_real(w_in)) or (not check_pos_real(bin_width)):
            raise ValueError("Input must be two positive real numbers")

    nw_pts = number_of_window_points(w_in, bin_width)
    return numpy.zeros(nw_pts) + 1.0
