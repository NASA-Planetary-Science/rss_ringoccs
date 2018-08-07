"""

calc_D_km.py

Purpose: Calculate the light distance between two input times.

Revisions:
    2018 Jul 09 - jfong - original
    2018 Aug 03 - jfong - remove for loop because inputs should be np.ndarrays,
                          not lists
"""

import numpy as np
import spiceypy as spice


def calc_D_km(t1, t2):
    """
    This calculates the light distance between two input times, in km.

    Args:
        t1 (np.ndarray):
            Array of time in seconds.
        t2 (np.ndarray): 
            Array of time in seconds.

    Returns:
        D_km_vals (np.ndarray): 
            Array of light distance in km.
    """
    D_km_vals = abs(t1-t2) * spice.clight()
    return D_km_vals
