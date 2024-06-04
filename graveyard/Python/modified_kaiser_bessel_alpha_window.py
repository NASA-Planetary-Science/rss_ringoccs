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
#       Computes the modified Kaiser Bessel window for arbitrary alpha.        #
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

def kbmdal(w_in, bin_width, alpha, check_error = True):
    """
        Function:
            kbmdal
        Purpose:
            Create the modified Kaiser-Bessel window
            with arbitrary alpha value.
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
                The modified Kaiser-Bessel alpha window of
                width w_in and spacing bin_width between points.
        Dependencies:
            [1] numpy
            [2] scipy.special
        Notes:
            [1] The Modified Kaiser-Bessel window is computed using
                the modified Bessel Function of the First Kind. It's
                value is:
                        I0(alpha * sqrt(1 - 4x^2 / w^2)) - 1
                    y = ------------------------------------
                                    I0(alpha) - 1
                where w is the window width.
            [2] The endpoints of the Modified Kaiser-Bessel function
                are zero. This means that the modified version has
                no discontinuities in it.
            [3] The Kaiser-Bessel functions and the modified
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
    factor = 1.0 / (scipy.special.iv(0.0, alpha) - 1.0)

    # The argument for the modified Bessel function.
    arg = alpha * numpy.sqrt(1.0 - numpy.square(2.0 * arr / w_in))

    # The window function.
    return (scipy.special.iv(0.0, arg) - 1.0) * factor
