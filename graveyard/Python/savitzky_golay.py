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
#       Smooths a function using the Savitzky-Golay filter.                    #
#   Notes:                                                                     #
#       This is a modification of the scipy routine.                           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import math
import numpy

def savitzky_golay(y_vals, window_size, order, deriv = 0):
    """
        Purpose:
            To smooth data with a Savitzky-Golay filter. This removes
            high frequency noise while maintaining many of the
            original features of the input data.
        Arguments:
            y_vals (numpy.ndarray):
                The input "Noisy" data.
            window_size (int):
                The length of the window. Must be an odd number.
            order (int):
                The order of the polynomial used for filtering.
                Must be less then window_size - 1.
        Keywords:
            deriv (int):
                The order of the derivative what will be computed.
        Output:
            y_smooth (numpy.ndarray):
                The data smoothed by Savitzky-Golay filter.
        Notes:
            The Savitzky-Golay is a type of low-pass filter, particularly
            suited for smoothing noisy data. The main idea behind this
            approach is to make for each point a least-square fit with a
            polynomial of high order over a odd-sized window centered at
            the point.
        References:
            [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
                Data by Simplified Least Squares Procedures. Analytical
                Chemistry, 1964, 36 (8), pp 1627-1639.
            [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
                W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
                Cambridge University Press ISBN-13: 9780521880688
    """

    # The previous routine when something along the lines as follows:
    #   order_range = numpy.arange(order + 1)
    #   half_window = (window_size -1) // 2
    # precompute coefficients
    #   b = numpy.mat(
    #       [
    #           [k**i for i in order_range]
    #           for k in numpy.arange(-half_window, half_window+1)
    #       ]
    #   )
    #   m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with values taken from the signal itself.
    #   firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    #   lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    #   y = numpy.concatenate((firstvals, y, lastvals))
    #   return numpy.convolve(m[::-1], y, mode='valid')

    if window_size < 1:
        raise ValueError("Window size must be positive.")

    if order < 1:
        raise ValueError("The polynomial order must be positive.")

    if window_size % 2 != 1:
        raise ValueError("Window size must be an odd integer.")

    if window_size < order + 2:
        raise ValueError("Window size must be less than the polynomial order.")

    half_window = (window_size - 1) // 2

    # precompute coefficients
    b_matrix = numpy.zeros((window_size, order+1))
    b_matrix[..., 0] = 1

    for k in range(half_window):
        current_index = (half_window) - k
        right_index = current_index
        left_index = -current_index

        for j in range(1, order+1):
            b_matrix[k, j] = left_index
            b_matrix[window_size-1-k, j] = right_index
            left_index *= -current_index
            right_index *= current_index

    b_matrix = numpy.mat(b_matrix)
    m_factor = math.factorial(deriv) * numpy.linalg.pinv(b_matrix).A[deriv]

    # Pad the endpoints with values from the signal.
    firstvals = y_vals[0]-numpy.abs(y_vals[1:half_window+1][::-1]-y_vals[0])
    lastvals = y_vals[-1]+numpy.abs(y_vals[-half_window-1:-1][::-1]-y_vals[-1])
    y_smooth = numpy.concatenate((firstvals, y_vals, lastvals))
    return numpy.convolve(m_factor[::-1], y_smooth, mode = 'valid')
