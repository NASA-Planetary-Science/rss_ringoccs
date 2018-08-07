"""
get_pole.py

Purpose: Get unit vector in pole direction from kernel constants.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def get_pole(et, planet, kernels=None):
    """
    This calculates unit vector in pole direction from kernel constants.

    Args:
        et (float64): 
            Ephemeris seconds past J2000
        planet (str):
            Planet name -- must be compatible with NAIF
        kernels (list): 
            List of NAIF kernels, including path

    Output:
        nhat_p (np.ndarray): 
            1x3 unit vector in pole direction.

    Note:
        [1] Quadratic terms for pole direction are typically zero but
            are retained here for consistency with PCK file format definitions.
    """

    # Load kernels
    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Retrieve right ascension and declination from kernel pool
    bodynm = planet
    item = 'POLE_RA'
    maxn = 3
    dim1, pole_RA = spice.bodvrd(bodynm, item, maxn)

    bodynm = planet
    item = 'POLE_DEC'
    maxn = 3
    dim2, pole_DEC = spice.bodvrd(bodynm, item, maxn)

    # Calculate pole ra and dec using quadratic terms
    dt_centuries = et / (spice.spd()*365.25*100.)
    rap = (pole_RA[0] + dt_centuries*pole_RA[1]
            + dt_centuries**2*pole_RA[2])
    dep = (pole_DEC[0] + dt_centuries*pole_DEC[1]
            + dt_centuries**2*pole_DEC[2])

    # Convert to rectangular coordinates
    inrange = 1.
    re = rap * spice.rpd()
    dec = dep * spice.rpd()
    nhat_p = spice.radrec(inrange, re, dec)

    return nhat_p
