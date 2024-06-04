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
#       Computes the normalized equivalent width of a window function.         #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""
import numpy
from .check_real import check_real

def compute_norm_eq(w_func):
    """
        Function:
            compute_norm_eq
        Purpose:
            Computes normalized equivalenth width of a given function.
        Variables:
            w_func:
                Any function (usually a window function).
        Outputs:
            normeq:
                The normalized equivalent width of w_func.
        Dependencies:
            [1] numpy
        Notes:
            The normalized equivalent width is effectively computed
            using Riemann sums to approximate integrals. Therefore
            large dx values (Spacing between points in w_func)
            will result in an inaccurate normeq. One should keep
            this in mind during calculations.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        History:
            2018/05/16: Ryan Maguire
                Created.
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
        """

    if not check_real(w_func):
        raise TypeError(w_func)

    nw_pts = float(numpy.size(w_func))
    numer = numpy.sum(numpy.square(w_func))
    denom = numpy.square(numpy.sum(w_func))
    return nw_pts * numer / denom
