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
#       Provides a class for working with window functions.                    #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
from .check_bool import check_bool
from .check_pos_real import check_pos_real

class WindowRangeData:
    """
        Class:
            WindowRangeData
        Purpose:
            Stores the data needed to computed the window range.
    """
    def __init__(self, normeq, fsky, fres, rho_dot):
        self.set_resolution(1.0)
        self.set_allen_deviation(2.0E-13)
        self.set_b_factor(True)
        self.set_normalized_equivalent_width(normeq)
        self.fsky = fsky
        self.fres = fres
        self.rho_dot = rho_dot

    def set_normalized_equivalent_width(self, normeq):
        """
            Sets the normalized equivalent width and checks for errors.
        """
        if not check_pos_real(normeq):
            raise ValueError("The equivalent width must be a positive number.")

        self.normeq = float(normeq)

    def set_resolution(self, res):
        """
            Sets the desired resolution and checks for errors.
        """
        if not check_pos_real(res):
            raise ValueError("The resolution must be a positive number.")

        self.res = float(res)

    def set_allen_deviation(self, sigma):
        """
            Sets the Allen deviation and checks for errors.
        """
        if not check_pos_real(sigma):
            raise ValueError("The Allen deviation must be a positive number.")

        self.sigma = float(sigma)

    def set_b_factor(self, b_factor):
        """
            Sets the b-factor for the inverse resolution function.
        """
        if not check_bool(b_factor):
            raise TypeError("b-factor must be a Boolean.")

        self.bfac = bool(b_factor)
