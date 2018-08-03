"""

calc_D_km.py

Purpose: Calculate the light distance between two input times.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_D_km(t1, t2):
    """
    This calculates the light distance between two input times, in km.

    Args:
        t1 (np.ndarray): Array of time in seconds.
        t2 (np.ndarray): Array of time in seconds.

    Returns:
        D_km_vals (np.ndarray): Array of light distance in km.
    """

    npts = len(t1)

    D_km_vals = np.zeros(npts)

    for n in range(npts):
        dt = abs(t1[n] - t2[n])
        D_km = dt * spice.clight()
        D_km_vals[n] = D_km

    return D_km_vals
