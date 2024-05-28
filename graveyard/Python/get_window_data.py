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
#       Computes the number of points in a window and creates an array.        #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
from .check_pos_real import check_pos_real
from .get_window_array import get_window_array
from .number_of_window_points import number_of_window_points

def get_window_data(w_in, bin_width, check_error = True):
    """
        Function:
            get_window_data
        Purpose:
            Computes variables that are common to many window
            functions. This includes the number of points in
            the window, and an array for the independent variable
            for the windowing function.
        Arguments:
            w_in (float):
                The width of the window. Should be a
                positive real number.
            bin_width (float):
                The distance between samples, or bins. Should
                also be a positive real number.
        Keywords:
            check_error (bool):
                Boolean for whether or not to check the inputs
                for errors before computing.
        Outputs:
            nw_pts (int):
                The number of points in the window.
            arr (numpy.ndarray):
                Array for the independent variable of the window.
        Called Functions:
            check_pos_real:
                Determines if a variable is a positive real number.
            get_window_array:
                Creates an array for the independent variable.
            number_of_window_points:
                Computes the number of points in the window.
    """
    if check_error:
        if (not check_pos_real(w_in)) or (not check_pos_real(bin_width)):
            raise ValueError("Input must be two positive real numbers")

    nw_pts = number_of_window_points(w_in, bin_width)
    arr = get_window_array(nw_pts, bin_width)

    return nw_pts, arr
