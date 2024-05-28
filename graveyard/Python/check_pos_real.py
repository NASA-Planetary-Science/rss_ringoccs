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
#       Checks if a variable is a positive real number.                        #
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

def check_pos_real(pos):
    """
        Function:
            check_pos_real
        Purpose:
            Check if an input is a single positive number.
        Arguments:
            pos:
                Any valid input (Int, float, string, list, etc).
        Outputs:
            pout (bool):
                A Boolean (True or False), depending on whether
                or not pos is a single positive real number.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [1] or
            [1.0], which looks positive and real, the function will
            return False.
        References:
            [1] John Stillwell, The Real Numbers:
                An Introduction to Set Theory and Analysis,
                Springer International Publishing, 2013
            [2] https://en.wikipedia.org/wiki/Real_number
        History:
            2018/05/14: Ryan Maguire
                Translated from IDL.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    if not type(pos) in VALID_TYPES.rtypes:
        return False

    if numpy.size(pos) != 1:
        return False

    if pos <= 0:
        return False

    return True
