"""

calc_sc_state.py

Purpose: Calculate spacecraft state vector in a planetocentric frame.

Revisions:
    2018 Jul 10 - jfong - original
"""

import numpy as np
import spiceypy as spice
from .xform_j2k_to_pcf import xform_j2k_to_pcf


def calc_sc_state(et_vals, spacecraft, planet, dsn, nhat_p, kernels=None):
    """
    This calculates spacecraft state vector in a planetocentric frame.

    Args:
        et_vals (np.ndarray): Array of spacecraft event times in ET seconds.
        spacecraft (str): Spacecraft name -- must be compatible with NAIF.
        planet (str): Planet name -- must be compatible with NAIF.
        dsn (str): DSN observing station ID -- must be compatible with NAIF.
        nhat_p (np.ndarray): 1x3 array unit vector in planet pole direction.
        kernels (list): List of NAIF kernels, including path.

    Output:
        R_sc_km_vals (list): List of np.ndarrays of spacecraft position vector
                        in km.
        R_sc_dot_kms_vals (list): List of np.ndarrays of spacecraft velocity
                        vector in km/s.
    """

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernel)

    npts = len(et_vals)
    R_sc_km_vals = []
    R_sc_dot_kms_vals = []
    for n in range(npts):
        et = et_vals[n]

        # Compute planet to spacecraft state vector in J2000 frame,
        #   with no light-time correction
        targ = spacecraft
        ref = 'J2000'
        abcorr = 'NONE'
        obs = planet
        starg0, ltime0 = spice.spkezr(targ, et, ref, abcorr, obs)

        R_sat2sc_km = starg0[0:3]
        R_sat2sc_dot_kms = starg0[3:6]

        # Transform vectors to planetocentric frame
        R_sc_km = xform_j2k_to_pcf(R_sat2sc_km, et, spacecraft, dsn, nhat_p)
        R_sc_dot_kms = xform_j2k_to_pcf(R_sat2sc_dot_kms, et,
                spacecraft, dsn, nhat_p)

        R_sc_km_vals.append(R_sc_km)
        R_sc_dot_kms_vals.append(R_sc_dot_kms)

    return R_sc_km_vals, R_sc_dot_kms_vals
