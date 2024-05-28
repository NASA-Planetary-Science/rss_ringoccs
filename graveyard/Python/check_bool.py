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
#       Checks if a variable is a Boolean.                                     #
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

def check_bool(boo):
    """
        Function:
            check_boole
        Purpose:
            Check if an input is Boolean.
        Variables:
            boo:
                Any valid input (Int, float, string, list, etc).
        Outputs:
            bout (bool):
                A Boolean (True or False), depending on whether
                or not boo is Boolean.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [0] or
            [1], which looks Boolean, the function will return False.
        References:
            [1] https://en.wikipedia.org/wiki/Boolean_data_type
        History:
            2018/05/14: Ryan Maguire
                Translated from IDL.
            2018/05/16: Ryan Maguire
                Allow True and False inputs (Boolean literals).
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    boo_type = type(boo)
    if not (boo_type in VALID_TYPES.rtypes) and (boo_type != type(True)):
        return False

    if numpy.size(boo) != 1:
        return False

    if (boo != 0) and (boo != 1) and (boo is not True) and (boo is not False):
        return False

    return True
