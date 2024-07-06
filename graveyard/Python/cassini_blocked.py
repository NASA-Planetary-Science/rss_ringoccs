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
#       Using spice calls to determine if Cassini is behind (1) Saturn or      #
#       Saturn's atmopshere, or (2) behind Saturn or Saturn's ionosphere at a  #
#       given set of times. Uses spiceypy.occult(), with Cassini as target and #
#       DSN station as observer.                                               #
################################################################################
#   Author:     Glenn Steranka                                                 #
#   Date:       2018/04/27                                                     #
################################################################################
"""
import numpy
import spiceypy
from . spm_to_et import spm_to_et


def cassini_blocked(spm_vals, rsr_inst, kernels, verbose = False):
    """
        Checks if Cassini is blocked by Saturn.
    """

    spiceypy.kclear()
    spiceypy.furnsh(kernels)

    et_vals = spm_to_et(spm_vals, rsr_inst.doy, rsr_inst.year, kernels=kernels)
    n_pts = len(et_vals)

    # Test to verify that we converted to ephemeris time correctly
    if verbose:
        print('First UTC time: ', spiceypy.et2utc(et_vals[0], 'D', 3))

    # CN aberration correction
    abcorr = 'CN'

    # Get original radii
    _, body699_radii = spiceypy.bodvrd('SATURN', 'RADII', 3)

    if verbose:
        print(
            'Original radii for Saturn: ' + str(body699_radii[0]) + ", " +
            str(body699_radii[1]) + ", " + str(body699_radii[2])
        )

    # Radii to look for atmosphere interference and ionosphere interference
    body699_atm = [rho699 + 500.0 for rho699 in body699_radii]
    body699_ion = [rho699 + 5000.0 for rho699 in body699_radii]

    # Use atmospheric radii in kernel pool
    spiceypy.pdpool('BODY699_RADII', body699_atm)

    if verbose:
        _, body699_atm = spiceypy.bodvrd('Saturn', 'RADII', 3)
        print(
            "Atmospheric radii: " + str(body699_atm[0]) + ", " +
            str(body699_atm[1]) + ", " + str(body699_atm[2])
        )

    # Look for atmospheric occultation
    occ_code_vals_atm = numpy.zeros(n_pts)
    for i in range(n_pts):
        occ_code_vals_atm[i] = spiceypy.occult(
            'Saturn', 'ELLIPSOID', 'IAU_SATURN', 'Cassini', 'POINT', ' ',
            abcorr, rsr_inst.dsn, et_vals[i]
        )

    is_blocked_atm = occ_code_vals_atm != 0

    # Use ionosphere radii in kernel pool
    spiceypy.pdpool('BODY699_RADII', body699_ion)

    if verbose:
        _, body699_ion = spiceypy.bodvrd('Saturn', 'RADII', 3)
        print(
            "Ionospheric radii: " + str(body699_ion[0]) + ", " +
            str(body699_ion[1]) + ", " + str(body699_ion[2])
        )

    # Look for ionospheric occultation
    occ_code_vals_ion = numpy.zeros(n_pts)

    for i in range(n_pts):
        occ_code_vals_ion[i] = spiceypy.occult(
            'Saturn', 'ELLIPSOID', 'IAU_SATURN', 'Cassini', 'POINT', ' ',
            abcorr, rsr_inst.dsn, et_vals[i]
        )

    is_blocked_ion = occ_code_vals_ion != 0

    if verbose:

        if numpy.any(is_blocked_atm):
            print(
                'First SPM occulted by atmosphere: ',
                (spm_vals[is_blocked_atm])[0]
            )

        else:
            print('No atmospheric occultation detected')

        if numpy.any(is_blocked_ion):
            print(
                'First SPM occulted by ionosphere: ',
                (spm_vals[is_blocked_ion])[0]
            )
        else:
            print('No ionospheric occultation detected')

    # Put regular Saturn radius into kernel pool
    spiceypy.pdpool('BODY699_RADII', body699_radii)

    # Test that we successfuly went back to original radii
    if verbose:
        _, body699_radii = spiceypy.bodvrd('Saturn', 'RADII', 3)
        print(
            'Returning to original Saturn radii: ' +
            str(body699_radii[0]) + ", " +
            str(body699_radii[1]) + ", " +
            str(body699_radii[2])
        )

    return is_blocked_atm, is_blocked_ion
