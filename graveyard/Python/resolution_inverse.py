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
#       Computes the inverse resolution function.                              #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
import numpy
import scipy.special
from .check_real import check_real
from .check_complex import check_complex


def resolution_inverse(x_val):
    """
        Function:
            resolution_inverse
        Purpose:
            Computes the inverse of y = x / (exp(-x) + x - 1).
        Variables:
            x_val (numpy.ndarray):
                A real or complex array.
        Outputs:
            out_val (numpy.ndarray):
                The inverse of x / (exp(-x) + x - 1)
        Dependencies:
            [1] numpy
            [2] scipy.special
        Notes:
            If the input is not in the atypes list, this function
            will raise a TypeError, returning the user to
            whichever shell they ran this function from. This
            function wants ints, floats, complex numbers, or arrays.
        Method:
            The inverse of x/(exp(-x)+x-1) is computed using the
            LambertW function. This function is the inverse of
            y = x * exp(x). This is computed using the scipy.special
            lambertw function.
        Warnings:
            [1] The real part of the argument must be greater than 1.
            [2] The scipy.special lambertw function is slightly
                inaccurate when it's argument is near -1/e. This
                argument is z = exp(x/(1-x)) * x/(1-x)
        References:
            [1] http://mathworld.wolfram.com/LambertW-Function.html
            [2] https://en.wikipedia.org/wiki/Lambert_W_function
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """

    # Booleans for the type of the input.
    is_real = check_real(x_val)
    is_complex = check_complex(x_val)

    # The function is only defined for real and complex variables.
    if not is_real and not is_complex:
        raise TypeError("Input must be real or complex.")

    # Intermediate factors to save on divisions.
    one_minus_x = 1.0 - x_val
    rcpr_one_minus_x = 1.0 / one_minus_x
    factor = x_val * rcpr_one_minus_x

    # The Lambert W factors for the inverse calculation.
    lambert_factor = numpy.exp(factor) * factor
    lambert_eval = scipy.special.lambertw(lambert_factor)

    # The formula for the inverse function.
    out_val = (one_minus_x * lambert_eval - x_val) * rcpr_one_minus_x

    # Real values map to real values. Return the real part if the input
    # has zero for the imaginary part.
    if is_real:
        return numpy.real(out_val)

    return out_val
