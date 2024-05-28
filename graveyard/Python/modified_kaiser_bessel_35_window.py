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
#       Computes the modified Kaiser Bessel window with alpha = 3.5 pi.        #
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

def kbmd35(w_in, bin_width, check_error = True):
    """
        Function:
            kbmd35
        Purpose:
            Create the Modified Kaiser-Bessel window with alpha = 3.5 pi.
        Arguments:
            w_in (float):
                Window width.
            bin_width (float):
                Width of one point (or one bin).
        Outputs:
            w_func (numpy.ndarray):
                The modified Kaiser-Bessel window with alpha = 3.5 pi.
        Dependencies:
            [1] numpy
            [2] scipy
        Notes:
            [1] The Modified Kaiser-Bessel window is computed using
                the modified Bessel Function of the First Kind. It's
                value is:
                        I0(alpha * sqrt(1 - 4x^2 / w^2)) - 1
                    y = ------------------------------------
                                    I0(alpha) - 1
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kbmd35 window function has an alpha value of
                alpha = 3.5 * pi
            [3] The endpoints of the Modified Kaiser-Bessel function
                are zero. This means that the modified version has
                no discontinuities in it.
            [4] The Kaiser-Bessel functions and the modified
                Kaiser-Bessel functions are equal at the center of
                the window.
        Warnings:
            [1] Unlike the Kaiser-Bessel functions, the modified
                Kaiser-Bessel functions evaluate to zero at the
                endpoints of the window.
            [2] Small alpha values will result in the Kaiser-Bessel
                function and the modified Kaiser-Bessel disagreeing
                dramatically. For example, alpha = 0 gives a constant
                curve for the Kaiser-Bessel, but a bell-shaped curve
                for the modified version.
        References:
            [1] https://en.wikipedia.org/wiki/Window_function
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

    # 1 / (I0(3.5 pi) - 1) where I0 is the zeroth modified Bessel function.
    factor = 1.3780177090677150E-04

    # The argument for the modified Bessel function.
    arg = alpha * numpy.sqrt(1.0 - numpy.square(2.0 * arr / w_in))

    # The window function.
    return (scipy.special.iv(0.0, arg) - 1.0) * factor
