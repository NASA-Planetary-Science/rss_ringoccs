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
#       Gets the current date as a string.                                     #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/27                                                         #
################################################################################
"""

import time

def date_string():
    """
        Function:
            date_string
        Purpose:
            Creates a date string of the form "yyyymmddHHMMSS"
        Arguments:
            None.
        Outputs:
            date (str):
                The current date "year/month/day/hour/minute/second"
        Dependencies:
            [1] time
    """

    # The time module has just the function for this. Simply need the format.
    return time.strftime("%Y%m%d%H%M%S")

def date_string_with_underscores():
    """
        Function:
            date_string_with_underscores
        Purpose:
            Creates a date string of the form "yyyy_mm_dd_HH_MM_SS"
        Arguments:
            None.
        Outputs:
            date (str):
                The current date "year/month/day/hour/minute/second"
        Dependencies:
            [1] time
    """
    return time.strftime("%Y_%m_%d_%H_%M_%S")
