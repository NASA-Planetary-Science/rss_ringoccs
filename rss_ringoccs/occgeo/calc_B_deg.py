"""

calc_B_deg.py

Purpose: Calculate ring opening angle, or the observed ring elevation.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_B_deg(et_vals, spacecraft, dsn, nhat_p, kernels=None):
    """
    This calculates ring opening angle, or the observed ring elevation.

    Args:
        et_vals (np.ndarray): Array of earth-received times in ET seconds.
        dsn (str): DSN observing station ID -- must be compatible with NAIF.
        nhat_p (np.ndarray): 1x3 array unit vector in planet pole direction.
        kernels (list): List of NAIF kernels, including path.

    Returns:
        B_deg_vals (np.ndarray): Array of ring opening angle in degrees.
    """

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernel)
    npts = len(et_vals)
    B_deg_vals = np.zeros(npts)

    for n in range(npts):
        # Compute Cassini to dsn position vector
        targ = dsn
        et = et_vals[n]
        ref = 'J2000'
        abcorr = 'CN'
        obs = spacecraft
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        # Calculate B as the complement to the angle made by the
        #   Saturn pole vector and the Cassini to DSN vector
        v1 = starg
        v2 = nhat_p
        B_rad = (np.pi/2.) - spice.vsep(v1, v2)

        B_deg_vals[n] = B_rad * spice.dpr()

    return B_deg_vals
