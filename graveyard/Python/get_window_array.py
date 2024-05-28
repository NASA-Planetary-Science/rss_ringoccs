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
#       Creates an array for the independent variable of a window function.    #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy:                                                                 #
#           Used for the arange function.                                      #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""
import numpy

def get_window_array(nw_pts, bin_width):
    """
        Function:
            get_window_array
        Purpose:
            Creates an array for the window from the number of
            points and displacement (bin width) between samples.
        Arguments:
            nw_pts (int):
                The number of points in the window.
            bin_width (float):
                Width of one point (or one bin).
        Output:
            arr (numpy.ndarray):
                An array for the window function.
        Notes:
            No error checks are performed. The inputs must be a
            positive integer and a positive floating point number.
    """
    return (numpy.arange(nw_pts) - (0.5 * (nw_pts - 1))) * bin_width
