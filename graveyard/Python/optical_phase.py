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
#       Computes the optical phase from a complex transmittance.               #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import numpy

def optical_phase(transmittance):
    """
        Function:
            optical_phase
        Purpose:
            Compute the phase from the complex transmittance.
        Variables:
            transmittance:
                The complex transmittance.
        Output:
            phase:
                The phase (in radians).
        Dependencies:
            [1] numpy
        Notes:
            [1] The phase of the complex transmittance is the angle
                made it makes with the x-axis in the complex plane.
                This is the arctangent of the ratio of the imaginary
                part to the real part. This equation works for both
                diffracted and reconstructed transmittances.
        Refernces:
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

    return numpy.angle(transmittance)
