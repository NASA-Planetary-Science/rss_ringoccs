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
#       Computes the Kaiser Bessel window function with real alpha value.      #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
#   2.) scipy                                                                  #
#           The scipy.special module contains Bessel functions that are used.  #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""
import numpy
import scipy.special
from .get_window_data import get_window_data

def kbal(w_in, bin_width, alpha, check_error = True):
    """
        Function:
            kbal
        Purpose:
            Create Kaiser-Bessel window with arbitrary alpha value.
        Variables:
            w_in (float):
                Window width.
            bin_width (float):
                The width of one bin.
            alpha (float):
                The alpha parameter for the Kaiser-Bessel window.
        Keywords:
            check_error (bool):
                Boolean for checking the inputs for errors.
        Outputs:
            w_func (numpy.ndarray):
                The Kaiser-Bessel alpha window of width w_in
                and spacing bin_width between points.
        Dependencies:
            [1] numpy
            [2] scipy.special
        Notes:
            [1] The Kaiser-Bessel window is computed using the
                modified Bessel Function of the First Kind. It's
                value is:
                    y = I0(alpha * sqrt(1 - 4x^2/w^2)) / I0(alpha)
                where w is the window width.
            [2] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1 + 2*alpha) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1/I_0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha values beyond
                2.5 * pi this effect is negligible.
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
            2018/05/16: Ryan Maguire
                Translated from IDL.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    _, arr = get_window_data(w_in, bin_width, check_error = check_error)

    # 1 / I0(alpha) where I0 is the zeroth modified Bessel function.
    factor = 1.0 / scipy.special.iv(0.0, alpha)

    # The argument for the modified Bessel function.
    arg = alpha * numpy.sqrt(1.0 - numpy.square(2.0 * arr / w_in))

    # The window function.
    return scipy.special.iv(0.0, arg) * factor
