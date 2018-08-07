"""

calc_elevation_deg.py

Purpose: Function to calculate the elevation of a target above the horizon for
         a given observer.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_elevation_deg(et_vals, target, obs, kernels=None):
    """
    Calculate the elevation of a target above the horizon for a given observer.

    Args:
        et_vals (np.ndarray): 
            Array of observed event times in ET sec.
        target (str):         
            Target name -- must be compatible with NAIF. This will typically
            be spacecraft or planet name.
        obs (str):
            Observer name -- must be compatible with NAIF. This will typically
            be observing dsn station name. Observer is assumed to be on Earth.
        kernels (list):       
            List of NAIF kernels, including path.

    Output:
        elev_deg_vals (np.ndarray): 
            Array of elevation angles in degrees.
    """

    npts = len(et_vals)
    elev_deg_vals = np.zeros(npts)

    # Load kernels
    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    for n in range(npts):
        et = et_vals[n]

        # Compute observer to target position vector in J2000
        #   with light-correction
        ref = 'J2000'
        abcorr = 'CN'
        ptarg1, ltime1 = spice.spkpos(target, et, ref, abcorr, obs)

        # Compute Earth to observer position vector in J2000
        #   without light correction
        abcorr = 'NONE'
        planet = 'EARTH'

        ptarg2, ltime2 = spice.spkpos(obs, et, ref, abcorr, planet)

        # Calculate elevation as the complement to the angle
        #   between ptarg1 (obs->target) and ptarg2 (Earth->obs)
        elev_deg_vals[n] = 90. - spice.vsep(ptarg1, ptarg2)*spice.dpr()

    return elev_deg_vals
