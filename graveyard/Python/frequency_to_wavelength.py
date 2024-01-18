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
#       Converts frequency (in Hertz) to wavelength (in kilometers).           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/14                                                         #
################################################################################
"""

def frequency_to_wavelength(frequency):
    """
        Function:
            frequency_to_wavelength
        Purpose:
            Converts frequency to wavelength, and vice-versa.
        Arguments:
            frequency:
                Frequency of the input in Hertz.
        Outputs:
            wavelength:
                Wavelength of the input in kilometers.
    """

    return 299792.0 / frequency
