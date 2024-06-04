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
#       Convert ephemeris time to seconds past midnight (SPM).                 #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/03/13                                                         #
################################################################################
"""
import spiceypy

def get_spm_from_et(et_val):
    """
        Function:
            get_spm_from_et
        Purpose:
            Converts a single ephemeris time to seconds past midnight.
    """

    # Convert ET to UTC string.
    utc_str = spiceypy.et2utc(et_val, "ISOD", 16)

    # Extract hour, minute, and second from UTC string.
    hour = float(utc_str[9:11])
    minute = float(utc_str[12:14])
    second = float(utc_str[15:])
    return hour*3600.0 + minute*60.0 + second
