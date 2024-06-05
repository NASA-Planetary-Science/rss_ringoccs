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
#       Gets a requested range and returns the equivalent list of values.      #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy
from .check_pos_real import check_pos_real
from .regions import REGIONS

def get_range_request(rngreq):
    """
        Function:
            get_range_request
        Purpose:
            This function takes a variety of input ranges and then
            converts them into the form [start,finish].
        Arguments:
            rngreq:
                The start/end points for diffraction correction.
                Preferred input is rngreq = [a, b]. Arrays are
                allowed and the range will be set as:
                    rng = [MIN(array),MAX(array)]
                Finally, certain strings containing a few of the
                regions of interests within the rings of Saturn
                are allowed. Strings are neither case nor space
                sensitive. For other planets use rng = [a, b]
        Output:
            rng:
                A list of the form [start,finish]
        Dependencies:
            [1] numpy
        Notes:
            [1] Inputs should be in kilometers. Numpy arrays, lists,
                and strings are allowed.
            [2] Strings are neither case nor space sensitive.
        References:
            [1] A.F. Cook, F.A. Franklin, F.D. Palluconi,
                Saturn's rings—A survey,
                Icarus, Volume 18, Issue 2, 1973, Pages 317-337,
                https://doi.org/10.1016/0019-1035(73)90214-5
            [2] Pollack, J.B., The Rings of Saturn,
                Space Science Reviews, Volume 18, Issue 1, 1975,
                Pages 3-93, https://doi.org/10.1007/BF00350197
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2024/06/05: Ryan Maguire
                Cleaned up and added to graveyard.
    """

    if isinstance(rngreq, list):
        if len(rngreq) != 2:
            raise TypeError("Requested range should be a list of two elements.")

        if not check_pos_real(rngreq[0]) or not check_pos_real(rngreq[1]):
            raise ValueError("Requested range should contain positive numbers.")

        # Make sure the entries are sorted correctly.
        if rngreq[0] < rngreq[1]:
            return rngreq

        # If not, flip the order.
        return [rngreq[1], rngreq[0]]

    # Several strings are allowed. These are provided by the REGIONS variable.
    if isinstance(rngreq, str):
        rng = rngreq.replace(" ", "").replace("'", "").replace("\"", "")
        rng = rng.lower()

        if rng not in REGIONS:
            raise ValueError("Invalid string for requested range.")

        return REGIONS[rng]

    if isinstance(rngreq, numpy.ndarray):
        minimum = numpy.min(rngreq)
        maximum = numpy.max(rngreq)

        if not check_pos_real(minimum) or not check_pos_real(maximum):
            raise ValueError("Requested range should contain positive numbers.")

        return [minimum, maximum]

    raise TypeError("Input should be a list or a string.")
