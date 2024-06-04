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
#       Provides a dictionary with different window functions.                 #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""

from .rectangular_window import rect
from .squared_cosine_window import coss
from .kaiser_bessel_20_window import kb20
from .kaiser_bessel_25_window import kb25
from .kaiser_bessel_35_window import kb35
from .modified_kaiser_bessel_20_window import kbmd20
from .modified_kaiser_bessel_25_window import kbmd25

# The dictionary consists of the window name, function, and normalized
# equivalent width.
FUNCTIONS = {
    "rect" : {"func" : rect, "normeq" : 1.00000000},
    "coss" : {"func" : coss, "normeq" : 1.50000000},
    "kb20" : {"func" : kb20, "normeq" : 1.49634231},
    "kb25" : {"func" : kb25, "normeq" : 1.65191895},
    "kb35" : {"func" : kb35, "normeq" : 1.92844639},
    "kbmd20" : {"func" : kbmd20, "normeq" : 1.52048174},
    "kbmd25" : {"func" : kbmd25, "normeq" : 1.65994218}
}
