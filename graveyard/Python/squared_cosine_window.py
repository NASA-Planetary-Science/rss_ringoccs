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
#       Provides the squared cosine window function.                           #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/15                                                         #
################################################################################
"""
import numpy
from .get_window_data import get_window_data

def coss(w_in, bin_width, check_error = True):
    """
        Function:
            coss
        Purpose:
            Create the cosine squared window.
        Arguments:
            w_in (float):
                Window width.
            bin_width (float):
                Width of one point (or one bin).
        Outputs:
            w_func (numpy.ndarray):
                The squared cosine window function of width w_in
                and spacing bin_width between points.
        Dependencies:
            [1] numpy
        Notes:
            [1] This window function is defined as:
                    y = cos^2(x * pi/w)
                where w is the window width.
            [2] The squared cosine window has the advantage of
                evaluating to zero at the endpoints of the window,
                meaning Gibb's effects from discontinuities do not
                arise like in the rectangular window and, to a much
                lesser degree, the Kaiser-Bessel windows.
            [3] The squared cosine window is wider than the all of
                the Kaiser-Bessel windows for alpha > 2.5*pi.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
        History:
            2018/05/15: Ryan Maguire
                Translated from IDL.
            2018/05/16: Ryan Maguire
                Made variables lower case.
            2024/05/28: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    _, arr = get_window_data(w_in, bin_width, check_error = check_error)
    return numpy.square(numpy.cos(numpy.pi * arr / w_in))
