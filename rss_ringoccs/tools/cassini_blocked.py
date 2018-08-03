"""

cassini_blocked.py

Purpose: Using spice calls to determine if Cassini is behind (1) Saturn or
         Saturn's atmopshere, or (2) behind Saturn or Saturn's ionosphere at a
         given set of times. Uses spice.occult(), with Cassini as target and
         DSN station as observer.

Revisions:
        gjs_cassini_behind_saturn_v2.py
    2018 Apr 27 - gsteranka - Original version
        gjs_cassini_behind_saturn_v4.py
    2018 May 03 - gsteranka - Original version
        cassini_blocked.py
    2018 May 03 - gsteranka - Original version
"""

import numpy as np
import pdb
import spiceypy as spice

try:
    from spm_to_et import spm_to_et
except ImportError:
    from .spm_to_et import spm_to_et


def cassini_blocked(spm_vals, rsr_inst, kernels, TEST=False):

    spice.kclear()
    spice.furnsh(kernels)

    et_vals = spm_to_et(spm_vals, rsr_inst.doy, rsr_inst.year, kernels=kernels)
    n_pts = len(et_vals)

    # Test to verify that we converted to ephemeris time correctly
    if TEST:
        UTCstart = spice.et2utc(et_vals[0], 'D', 3)
        print('First UTC time: ', UTCstart)

    # CN aberration correction
    abcorr = 'CN'

    # Get original radii
    dim, BODY699_RADII = spice.bodvrd('SATURN', 'RADII', 3)

    if TEST:
        print('Original radii for Saturn: %f, %f, %f' %
            (BODY699_RADII[0], BODY699_RADII[1], BODY699_RADII[2]))

    # Radii to look for atmosphere interference and ionosphere interference
    BODY699_ATM = [rho699 + 500.0 for rho699 in BODY699_RADII]
    BODY699_ION = [rho699 + 5000.0 for rho699 in BODY699_RADII]

    # Use atmospheric radii in kernel pool
    spice.pdpool('BODY699_RADII', BODY699_ATM)
    if TEST:
        dim, BODY699_ATM = spice.bodvrd('Saturn', 'RADII', 3)
        print('Atmospheric radii: %f, %f, %f' %
            (BODY699_ATM[0], BODY699_ATM[1],
            BODY699_ATM[2]))

    # Look for atmospheric occultation
    occ_code_vals_atm = np.zeros(n_pts)
    for i in range(n_pts):
        occ_code_vals_atm[i] = spice.occult('Saturn', 'ELLIPSOID', 'IAU_SATURN',
            'Cassini', 'POINT', ' ', abcorr, rsr_inst.dsn, et_vals[i])
    is_blocked_ATM = (occ_code_vals_atm != 0)

    # Use ionosphere radii in kernel pool
    spice.pdpool('BODY699_RADII', BODY699_ION)
    if TEST:
        dim, BODY699_ION = spice.bodvrd('Saturn', 'RADII', 3)
        print('Ionospheric radii: %f, %f, %f' %
            (BODY699_ION[0], BODY699_ION[1], BODY699_ION[2]))

    # Look for ionospheric occultation
    occ_code_vals_ion = np.zeros(n_pts)
    for i in range(n_pts):
        occ_code_vals_ion[i] = spice.occult('Saturn', 'ELLIPSOID', 'IAU_SATURN',
            'Cassini', 'POINT', ' ', abcorr, rsr_inst.dsn, et_vals[i])
    is_blocked_ION = (occ_code_vals_ion != 0)

    if TEST & np.any(is_blocked_ATM):
        print('First SPM occulted by atmosphere: ',
            (spm_vals[is_blocked_ATM])[0])
    elif TEST:
        print('No atmospheric occultation detected')

    if TEST & np.any(is_blocked_ION):
        print('First SPM occulted by ionosphere: ',
            (spm_vals[is_blocked_ION])[0])
    elif TEST:
        print('No ionospheric occultation detected')

    # Put regular Saturn radius into kernel pool
    spice.pdpool('BODY699_RADII', BODY699_RADII)

    # Test that we successfuly went back to original radii
    if TEST:
        dim, BODY699_RADII = spice.bodvrd('Saturn', 'RADII', 3)
        print('Returning to original Saturn radii: %f, %f, %f' %
            (BODY699_RADII[0], BODY699_RADII[1], BODY699_RADII[2]))

    return is_blocked_ATM, is_blocked_ION


if __name__ == '__main__':
    pass
