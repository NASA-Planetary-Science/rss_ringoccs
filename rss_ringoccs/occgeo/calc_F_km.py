"""

calc_F_km.py

Purpose: This calculates Fresnel scale using Eq. 6 of MTR1986.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_F_km(D_km_vals, f_sky_hz_vals, B_deg_vals, phi_ora_deg_vals):
    """
    This calculates the Fresnel scale using Eq. 6 of MTR1986.

    Args:
        D_km_vals (np.ndarray): Array of spacecraft to ring intercept point
                                distances in km.
        f_sky_hz_vals (np.ndarray): Array of downlink sinusoidal signal
                                frequency at the front-end of DSN station.
        B_deg_vals (np.ndarray): Array of ring opening angle in degrees.
        phi_ora_deg_vals (np.ndarray): Array of observed ring azimuth in
                                        degrees.
    Output:
        F_km_vals (np.ndarray): Array of Fresnel scale in km.

    Notes:
        [1] diffcorr uses an independently-calculated Fresnel scale

    """

    lambda_sky = spice.clight() / f_sky_hz_vals
    phi_ora_rad = spice.rpd() * phi_ora_deg_vals
    B_rad = spice.rpd() * B_deg_vals
    D_km = D_km_vals

    # Eq. 6 of MTR1986
    F_km_vals = np.sqrt((0.5 * lambda_sky * D_km * (1 - (np.cos(B_rad))**2
        * (np.sin(phi_ora_rad))**2)) / (np.sin(B_rad))**2)

    return F_km_vals
