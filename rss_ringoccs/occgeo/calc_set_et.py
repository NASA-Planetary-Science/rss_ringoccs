"""

calc_set_et.py

Purpose: Calculate spacecraft event time given observed event time.

Revisions:
    2018 Jul 10 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_set_et(et_vals, spacecraft, dsn, kernels=None):
    """
    This calculates the time at which photons left the spacecraft, in ET sec.

    Args:
        et_vals (np.ndarray): Array of earth-received times in ET seconds.
        spacecraft (str): Spacecraft name -- must be compatible with NAIF.
        dsn (str): DSN observing station ID -- must be compatible with NAIF

    Output:
        t_set_et_vals (np.ndarray): Array of spacecraft event times in ET sec.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Convert spacecraft and dsn strings to NAIF integer codes
    sc_code = spice.bodn2c(spacecraft)
    dsn_code = spice.bodn2c(dsn)

    # Use light travel time to calculate spacecraft event time
    t_set_et_vals = []
    for et in et_vals:
        set_et, ltime = spice.ltime(et, dsn_code, "<-", sc_code)
        t_set_et_vals.append(set_et)

    return np.asarray(t_set_et_vals)
