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
#       Provides a simple implementation of the Fresnel cosine function.       #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
#   2.) scipy                                                                  #
#           The scipy.special module contains the error function.              #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/14                                                         #
################################################################################
"""

# pylint incorrectly warns that scipy.special does not contain the erf
# function. Disable this warning.
# pylint: disable = no-member
import numpy
import scipy.special

def fresnel_cos(x_in):
    """
        Function:
            fresnel_cos
        Purpose:
            Compute the Fresnel cosine function.
        Variables:
            x_in:
                A real or complex argument.
        Outputs:
            f_cos:
                The Fresnel cosine integral of x_in.
        Dependences:
            [1] numpy
            [2] scipy.special
        Notes:
            [1] The Fresnel Cosine integral is the solution
                to the equation:
                    dy/dx = cos(pi/2 * x^2), y(0) = 0
                In other words:
                             -
                            | | x
                            |
                    C(x) =  | cos(pi t^2 / 2) dt
                          | |
                           - 0
                This is very closely related to the error function.
                Because of this, we may compute using:
                    C(x) = A erf(bx) + B erf(cx)
                where A, B, a, and b are complex constants.
            [2] The Fresnel Cosine and Sine integrals are computed by
                using the scipy.special Error Function. The Error
                Function, usually denoted Erf(x), is the solution to:
                    dy/dx = (2/sqrt(pi)) * exp(-x^2), y(0) = 0
                That is:
                                         -
                                        | | x
                                2       |
                    Erf(x) = --------   | exp(-t^2) dt
                             sqrt(pi) | |
                                       - 0
                Using Euler's Formula for exponentials allows one
                to use this to solve for the Fresnel Cosine integral.
            [3] The Fresnel Cosine integral is used for the solution
                of diffraction through a square well. Because of this
                it is useful for forward modeling problems in
                radiative transfer and diffraction.
        References:
            [1] https://en.wikipedia.org/wiki/Fresnel_integral
            [2] https://en.wikipedia.org/wiki/Error_function
            [3] http://mathworld.wolfram.com/FresnelIntegrals.html
            [4] http://mathworld.wolfram.com/Erf.html
        History:
            2018/05/14: Ryan Maguire
                Translated from IDL.
            2018/05/16: Ryan Maguire
                Allow complex arguments.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """

    # sqrt(pi) / 2 times the input.
    factor = 0.8862269254527579 * x_in

    # Multiplicative factors for the two error functions.
    left_scale = 0.25 - 0.25j
    right_scale = 0.25 + 0.25j

    # The Fresnel cosine decomposes as the sum of two error functions.
    left = left_scale * scipy.special.erf((1.0 + 1.0j) * factor)
    right = right_scale * scipy.special.erf((1.0 - 1.0j) * factor)
    f_cos = left + right

    # If the input is real, we'll return a real-valued output.
    is_real = numpy.isreal(x_in)

    # If the input is a single number, the .all() method is not defined.
    if numpy.size(is_real) == 1:
        if is_real:
            return numpy.real(f_cos)

        return f_cos

    # If the input is an array of real numbers, return a real-valued array.
    if is_real.all():
        return numpy.real(f_cos)
    return f_cos
