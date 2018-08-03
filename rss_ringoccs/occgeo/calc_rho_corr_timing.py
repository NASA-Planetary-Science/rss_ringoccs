"""

calc_rho_corr_timing.py

Purpose: Calculate radius correction due to timing offset with no slope term.

Revisions:
    2018 Jul 09 - jfong - original
"""


def calc_rho_corr_timing(rho_dot_kms, delta_t):
    """
    Calculate radius correction due to timing offset using no slope term.

    Args:
        rho_dot_kms (float64): Ring intercept radial velocity (array of value)
                                in km/s.
        delta_t (float64):     Time correction in seconds.

    Output:
        rho_corr_km (float64): Radius correction in km to be added to the
                               uncorrected radius.
    """
    rho_corr_km = rho_dot_kms * delta_t
    return rho_corr_km
