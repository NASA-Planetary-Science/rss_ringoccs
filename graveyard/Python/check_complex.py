
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
#       Checks if a variable is complex.                                       #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy:                                                                 #
#           Used for the size function.                                        #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/14                                                         #
################################################################################
"""
import numpy
from .valid_types import VALID_TYPES

def check_complex(com):
    """
        Function:
            check_complex
        Purpose:
            Check if an input is purely complex valued.
        Variables:
            com:
                Any valid input (Int, float, string, list, etc).
        Outputs:
            cout (bool):
                A Boolean (True or False), depending on whether
                or not com is complex valued.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the ctypes list, this function
            returns False. Even if the input is a list like [1j]
            or [1+1j], which looks complex valued, the function
            will return False.
        References:
            [1] Paul Nahin, An Imaginary Tale: The Story of i,
                Princeton Publishing, 1998
            [2] https://en.wikipedia.org/wiki/Complex_number
        History:
            2018/05/14: Ryan Maguire
                Translated from IDL.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    if not type(com) in VALID_TYPES.ctypes:
        return False

    if numpy.size(com) == 1:
        if numpy.imag(com) == 0:
            return False

        return True

    if (numpy.imag(com) == 0).all():
        return False

    return True
