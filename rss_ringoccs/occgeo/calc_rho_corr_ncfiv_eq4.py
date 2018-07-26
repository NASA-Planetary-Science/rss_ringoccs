"""

calc_rho_corr_ncfiv_eq4.py

Purpose: Calculate radius correction due to timing offset using Eq 4 of
         FRENCHETAL2017 NCF-IV paper.

Revisions:
    2018 Jul 09 - jfong - original
"""


def calc_rho_corr_ncfiv_eq4(rho_km, rho_dot_kms, delta_t, alpha):
    """
    Calculate radius correction due to timing offset using Eq4 of NCF-IV
        FRENCHETAL2017.

    Args:
        rho_km (float64):      Uncorrected radius (array of value) in km.
        rho_dot_kms (float64): Ring intercept radial velocity (array of value)
                                in km/s.
        delta_t (float64):     Time correction in seconds.
        alpha (float64):       Slope term.

    Output:
        rho_corr_km (float64): Radius correction in km to be added to the
                               uncorrected radius.
    """
    r0 = 100000.
    term1 = rho_dot_kms * delta_t
    term2 = alpha * ((rho_km - r0)/1000.)
    drho = term1 - term2
    return drho
