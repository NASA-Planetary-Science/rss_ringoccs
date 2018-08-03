"""

calc_beta.py

Purpose: Calculate B_eff_deg, the effective ring opening angle in degrees,
         and beta, the optical depth enhancement factor.

Revisions:
    2018 Jul 16 - jfong - original
"""

import numpy as np


def calc_B_eff_deg(B_deg, phi_ora_deg):
    """
    This calculates the effective ring opening angle in degrees.

    Args:
        B_deg (float64): Ring opening angle in degrees.
        phi_ora_deg (float64): Observed ring azimuth in degrees.

    Output:
        B_eff_deg (float64): Effective ring opening angle in degrees.

    References:
        [1] Gresh et al. (1986) Icarus 68, 481-502 Eq. 16 "An analysis of
            bending waves in Saturn's rings using Voyager radio occultation
            data"
    """
    B_rad = np.radians(B_deg)
    phi_ora_rad = np.radians(phi_ora_deg)
    B_eff_deg = np.arctan2(np.tan(B_rad), np.cos(phi_ora_rad))
    return B_eff_deg


def calc_beta(B_deg, phi_ora_deg):
    """
    This calculates the optical depth enhancement factor.

    Args:
        B_deg (float64): Ring opening angle in degrees.
        phi_ora_deg (float64): Observed ring azimuth in degrees.

    Output:
        beta (float64): Optical depth enhancement factor.
    """
    B_eff_deg = calc_B_eff_deg(B_deg, phi_ora_deg)
    beta = 1./(np.tan(B_eff_deg))
    return beta
