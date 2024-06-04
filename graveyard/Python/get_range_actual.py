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
#       Gets the allowed range of processing via window widths and data range. #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
import numpy
from .check_real import check_real
from .check_pos_real import check_pos_real

def get_range_actual(rho, rng, w_vals):
    """
        Function:
            get_range_actual
        Purpose:
            Compute the possible allowed range for processing, taking
            into consideration available data (rho) and the requested region.
        Variables:
            rho:
                Radial range of the data.
            rng:
                Requested start / end points for processing.
            w_vals:
                Window width as a function of ring radius.
        Output:
            start:
                The allowed starting point for processing.
            n_used:
                The number of points allowed for processing.
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
    """
    if not check_real(rho):
        raise TypeError("Rho must be an array of real numbers.")

    if numpy.min(rho) < 0.0:
        raise ValueError("Rho must be positive.")

    if numpy.size(rng) != 2:
        raise TypeError("Range must have format rng = [a, b].")

    if not check_pos_real(numpy.min(rng)):
        raise ValueError("Range must be positive.")

    if not check_real(w_vals):
        raise TypeError("w_vals must be real")

    if numpy.min(w_vals) < 0.0:
        raise ValueError("w_vals must be positive")

    if numpy.min(rng) > numpy.max(rho):
        raise ValueError("Requested range GREATER than available data.")

    if numpy.max(rng) < numpy.min(rho):
        raise ValueError("Requested range LESS than available data.")

    w_max = numpy.max(w_vals)
    rho_min_lim = numpy.min(rho)+numpy.ceil(w_max/2.0)
    rho_max_lim = numpy.max(rho)-numpy.ceil(w_max/2.0)
    rho_start = rho[numpy.min((rho >= numpy.min(rng)).nonzero())]
    rho_end = rho[numpy.max((rho <= numpy.max(rng)).nonzero())]
    rho_min = numpy.max([rho_min_lim,rho_start])
    rho_max = numpy.min([rho_max_lim,rho_end])
    start = int(numpy.min((rho >= rho_min).nonzero()))
    finish = int(numpy.max((rho <= rho_max).nonzero()))
    n_used = 1 + (finish - start)
    return start, n_used
