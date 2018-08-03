"""

calc_rho_km.py

Purpose: Calculate rho, the distance between Saturn center to ring intercept
         point, in km. This function returns only the scalar value of rho.
         To get the vector value, use function calc_rho_vec_km().

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice
from .calc_rho_vec_km import calc_rho_vec_km


def calc_rho_km(et_vals, planet, spacecraft, dsn, kernels=None, verbose=False):
    """
    Calculate the distance between Saturn center to ring intercept point.

    Args:
        et_vals (np.ndarray): Array of observed event times in ET sec.
        planet (str):         Planet name -- must be compatible with NAIF
        spacecraft (str):     Spacecraft name -- must be compatible with NAIF.
        dsn (str):            DSN observing station ID -- must be compatible
                              with NAIF.
        kernels (list):       List of NAIF kernels, including path.
        verbose (bool):       Default False prints nothing, while setting as
                              True prints steps

    Output:
        rho_km_vals (np.ndarray): Array of ring radii in km.
    """

    # Calculate rho vector
    rho_vec_km, t_ret_et = calc_rho_vec_km(et_vals, planet, spacecraft, dsn,
            kernels=kernels)

    rho_km_vals = []

    # Compute magnitude of each vector
    for vec in rho_vec_km:
        rho_km_vals.append(spice.vnorm(vec))

    return np.asarray(rho_km_vals)
