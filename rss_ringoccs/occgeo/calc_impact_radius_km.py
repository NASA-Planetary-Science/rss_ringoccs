"""

calc_impact_radius_km.py

Purpose: Calculate the radius of a sphere centered at Saturn and is tangent
         to the line-of-sight from spacecraft (at SET) to Earth (at OET).

Revisions:
    2018 Jul 10 - jfong - original
"""

import numpy as np
import spiceypy as spice
from .xform_j2k_to_pcf import xform_j2k_to_pcf


def calc_impact_radius_km(R_sc_km_vals, et_vals, spacecraft, dsn, nhat_p,
        kernels=None):
    """
    This calculates the impact radius in km.

    Args:
        R_sc_km_vals (list): List of 3-element arrays of spacecraft position
                vector in planetocentric frame at input et_vals.
        et_vals (np.ndarray): Array of Earth-received times in ET seconds.
        spacecraft (str): Spacecraft name -- must be compatible with NAIF.
        dsn (str): DSN observing station ID -- must be compatible with NAIF.
        nhat_p (np.ndarray): 1x3 array unit vector in planet pole direction.
        kernels (list): List of NAIF kernels, including path.

    Output:
        R_imp_km_vals (np.ndarray): Array of impact radius in km.
    """

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernel)

    npts = len(et_vals)
    R_imp_km_vals = np.zeros(npts)

    for n in range(npts):
        R_sc_km = R_sc_km_vals[n]
        et = et_vals[n]

        # Compute spacecraft to dsn position vector in J2000 frame,
        #   at et+ltime
        targ = dsn
        ref = 'J2000'
        abcorr = 'XCN'
        obs = spacecraft
        starg1, ltime1 = spice.spkpos(targ, et, ref, abcorr, obs)

        # Transform vector to planetocentric frame
        R_sc2dsn_km_pcf = xform_j2k_to_pcf(starg1, et, spacecraft,
                dsn, nhat_p)

        # Calculate distance from saturn center to point of closest approach
        lindir = R_sc2dsn_km_pcf
        linpt = R_sc_km
        point = [0., 0., 0.]
        pnear, distance = spice.nplnpt(linpt, lindir, point)

        R_imp_km_vals[n] = distance
    return R_imp_km_vals
