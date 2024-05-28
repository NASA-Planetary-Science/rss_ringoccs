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
#       Computes the number of samples in a window.                            #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy:                                                                 #
#           Used for the floor function.                                       #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/14                                                         #
################################################################################
"""
import numpy

def number_of_window_points(w_in, bin_width):
    """
        Function:
            number_of_window_points
        Purpose:
            Computes the number of points in a window.
        Arguments:
            w_in (float):
                Window width.
            bin_width (float):
                Width of one point (or one bin).
        Output:
            nw_pts (int):
                The number of points in the window.
        Notes:
            No error checks are performed. The inputs must be positive
            floating point numbers.
        History:
            2018/05/14: Ryan Maguire
                Original draft.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    return int(2.0 * numpy.floor(w_in / (2.0 * bin_width)) + 1.0)
