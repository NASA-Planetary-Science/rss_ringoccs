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
#       Provides a class with lists of various data types for real and complex #
#       numbers. These are used in error checking routines.                    #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/21                                                         #
################################################################################
"""
import numpy

class ValidTypes:
    """
        Class:
            ValidTypes
        Purpose:
            Contains listls of legal input types allowed in the
            various functions defined in rss_ringoccs.
        Arguments:
            There are no arguments to this class.
        Outputs:
            rtypes (list):
                A list of allowed data types that are real valued.
                Int, Float, and numpy arrays are included.
            ctypes (list):
                A list of allowed variable types that are complex valued.
                Complex and complex numpy arrays are included.
            atypes (list):
                List of all types. This is the concatenation of
                rtypes and ctypes.
        Notes:
            Lists are not included in any of the output lists.
            Some functions in rss_ringoccs allow lists as inputs,
            such as function pertaining to a requested range, but
            numpy arrays are also allowed.
            When in doubt, use numpy arrays.
        History:
            2018/05/21: Ryan Maguire
                Made into a class.
        """
    def __init__(self):
        """
            Method:
                __init__
            Purpose:
                Creates an instance of the class.
        """

        # Real types.
        self.rtypes = self.get_real_types()
        self.ctypes = self.get_complex_types()
        self.atypes = self.rtypes + self.ctypes

    def get_real_types(self):
        """
            Method:
                get_real_types
            Purpose:
                Returns a list of various real-valued data types.
            Arguments:
                None.
            Outputs:
                rtypes (list):
                    A list of real-valued data types.
        """
        return [
            type(0),
            type(0.0),
            type(numpy.intc(0)),
            type(numpy.intp(0)),
            type(numpy.int8(0)),
            type(numpy.int16(0)),
            type(numpy.int32(0)),
            type(numpy.int64(0)),
            type(numpy.array(0)),
            type(numpy.float16(0.0)),
            type(numpy.float32(0.0)),
            type(numpy.float64(0.0)),
            type(numpy.float128(0.0)),
            type(numpy.array(range(2)))
        ]

    def get_complex_types(self):
        """
            Method:
                get_complex_types
            Purpose:
                Returns a list of various complex-valued data types.
            Arguments:
                None.
            Outputs:
                ctypes (list):
                    A list of complex-valued data types.
        """
        return [
            type(complex(0)),
            type(numpy.complex64(0)),
            type(numpy.complex128(0)),
            type(numpy.complex256(0))
        ]

# Instance of the ValidTypes class.
VALID_TYPES = ValidTypes()
