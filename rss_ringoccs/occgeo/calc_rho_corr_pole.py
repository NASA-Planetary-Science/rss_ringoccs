"""

calc_rho_corr_pole.py

Purpose: Calculate radius correction due to an improved pole.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice
from .calc_rho_km import calc_rho_km


def calc_rho_corr_pole(rho_km_vals_ori, et_vals, dsn, planet, spacecraft,
        new_kernels):
    """
    Calculate radius correction due to improved pole.

    Args:
        rho_km_vals_ori (np.ndarray): Array of uncorrected ring radii in km.
        et_vals (np.ndarray): Array of earth-received times in ET sec
                              corresponding to rho_km_vals_ori.
        dsn (str): DSN observing station ID -- must be compatible with NAIF.
        planet (str): Planet name -- must be compatible with NAIF.
        spacecraft (str):     Spacecraft name -- must be compatible with NAIF.
        new_kernels (list): NAIF pck kernel with the improved pole.

    Output:
        rho_corr_km (np.ndarray): Radius corrections in km to be added to the
                               uncorrected radii.
    """

    # Calculate rho using new kernels
    rho_km_vals_new = calc_rho_km(et_vals, planet, spacecraft, dsn,
            kernels=new_kernels)

    # Radius correction is the difference between new and original radius
    rho_corr_km = rho_km_vals_new - rho_km_vals_ori
    return rho_corr_km
