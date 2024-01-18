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
#       Computes the optical depth from a complex transmittance.               #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import numpy

def optical_depth(transmittance, mu_vals):
    """
        Function:
            optical_depth
        Purpose:
            Compute the normalized optical depth.
        Variables:
            transmittance:
                The complex transmittance.
            mu_vals:
                The sine of the ring opening angle.
        Output:
            tau:
                The normalized optical depth.
        Dependencies:
            [1] numpy
        Notes:
            [1] The optical depth is proportional to the natural
                logarithm of the transmitted power. It is a unitless
                variable that arises in the study of radiative
                transfer, in particular the transfer equation. The
                normalized optical depth is normalized version of
                this, taking geometrical factors into account and
                using the normalized power.
            [2] Tau is often used to represent the optical depth. The
                equation used it tau = -mu * ln(Power), where mu is
                the sine of the ring opening angle (Denoted B), and
                where ln is the natural logarithm.
        References:
            [1] George B. Rybicki and Alan P. Lightman,
                Radiative Processes in Astrophysics,
                Wiley, 29 December 2007
            [2] S. Chandrasekhar, Radiative Transfer,
                Dover Publications, 1960
            [3] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [4] S. W. Asmar, R. G. French, E. A. Marouf, P. Schinder,
                J. W. Armstrong, P. Tortora, L. Iess, A. Anabtawi,
                A. J. Kliore, M. Parisi, M. Zannoni, and D. Kahan,
                Cassini Radio Science User's Guide, September, 2014,
                https://pds-rings.seti.org/cassini/rss/
    """

    # From log rules, ln(|P|^2) = 2 ln(|P|). We can save a call to numpy.square.
    return -2.0 * mu_vals * numpy.log(numpy.abs(transmittance))
