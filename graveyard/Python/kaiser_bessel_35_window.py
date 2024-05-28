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
#       Computes the Kaiser Bessel window function with alpha = 3.5 pi.        #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
#   2.) scipy                                                                  #
#           The scipy.special module contains Bessel functions that are used.  #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
import numpy
import scipy.special
from .get_window_data import get_window_data

def kb35(w_in, bin_width, check_error = True):
    """
        Function:
            kb35
        Purpose:
            Create the Kaiser-Bessel window with alpha = 3.5 pi.
        Arguments:
            w_in (float):
                Window width.
            bin_width (float):
                Width of one point (or one bin).
        Outputs:
            w_func (numpy.ndarray):
                The Kaiser-Bessel window with alpha = 3.5 pi.
        Dependencies:
            [1] numpy
            [2] scipy
        Notes:
            [1] The Kaiser-Bessel window is computed using the
                modified Bessel Function of the First Kind. It's
                value is:
                    y = I0(alpha * sqrt(1 - 4x^2/w^2)) / I0(alpha)
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kb35 window function has an alpha value of
                alpha = 3.5 * pi
            [3] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1 + 2*alpha) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1 / I0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha = 3.5 * pi,
                the endpoint is 0.00013778, which is small.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
            [3] On approximating the Modified Bessel Function of the
                First Kind and Toader-Qi Mean, Yang, ZH. & Chu, YM.
                J Inequal Appl (2016): 40., Springer,
                https://doi.org/10.1186/s13660-016-0988-1
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2018/05/16: Ryan Maguire
                Made variables lower case.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    _, arr = get_window_data(w_in, bin_width, check_error = check_error)

    # 3.5 pi, the parameter for the Kaiser-Bessel window.
    alpha = 10.995574287564276

    # 1 / I0(3.5 pi) where I0 is the zeroth modified Bessel function.
    factor = 1.3778278419510897E-04

    # The argument for the modified Bessel function.
    arg = alpha * numpy.sqrt(1.0 - numpy.square(2.0 * arr / w_in))

    # The window function.
    return scipy.special.iv(0.0, arg) * factor
