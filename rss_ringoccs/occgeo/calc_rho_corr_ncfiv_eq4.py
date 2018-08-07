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
        rho_km (float64 or np.ndarray):      
            Uncorrected radius in km.
        rho_dot_kms (float64 or np.ndarray): 
            Ring intercept radial velocity in km/s.
        delta_t (float64):
            Time correction in seconds.
        alpha (float64):       
            Slope term.

    Output:
        rho_corr_km (float64 or np.ndarray): 
            Radius correction in km to be added to the uncorrected radius.

    References:
        [1] Richard G. French, Colleen A. McGhee-French, Katherine Lonergan, 
            Talia Sepersky, Robert A. Jacobson, Philip D. Nicholson, 
            Mathew M. Hedman, Essam A. Marouf, Joshua E. Colwell,
            Noncircular features in Saturn’s rings IV: Absolute radius scale
            and Saturn’s pole direction, 
            Icarus, Volume 290, 2017, Pages 14-45, ISSN 0019-1035,
            https://doi.org/10.1016/j.icarus.2017.02.007.
        
    """
    r0 = 100000.
    term1 = rho_dot_kms * delta_t
    term2 = alpha * ((rho_km - r0)/1000.)
    drho = term1 - term2
    return drho
