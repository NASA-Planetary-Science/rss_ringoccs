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
#       Computes the optical power from a complex transmittance.               #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import numpy

def optical_power(transmittance):
    """
        Function:
            optical_power
        Purpose:
            Computes the power from complex transmittance.
        Arguments:
            transmittance:
                The complex transmittance.
        Output:
            power:
                The power.
        Dependencies:
            [1] numpy
        Notes:
            [1] The power is simply the square of the absolute value
                of the complex transmittance. This equation works for
                both a diffracted and a reconstructed transmittance.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] S. W. Asmar, R. G. French, E. A. Marouf, P. Schinder,
                J. W. Armstrong, P. Tortora, L. Iess, A. Anabtawi,
                A. J. Kliore, M. Parisi, M. Zannoni, and D. Kahan,
                Cassini Radio Science User's Guide, September, 2014,
                https://pds-rings.seti.org/cassini/rss/
    """

    return numpy.square(numpy.abs(transmittance))
