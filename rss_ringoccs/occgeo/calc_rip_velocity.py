"""

calc_rip_velocity.py

Purpose: Calculate ring intercept radial and azimuthal velocity.

Revisions:
    2018 Jul 10 - jfong - original
"""
import numpy as np


def calc_rip_velocity(rho_km_vals, phi_rl_deg_vals, dt):
    """
    This calculates the ring intercept point radial and azimuthal velocity.

    Args:
        rho_km_vals (np.ndarray): 
            Array of ring intercept points in km.
        phi_rl_deg_vals (np.ndarray): 
            Array of ring longitudes in degrees.
        dt (float64): 
            Constant time spacing between points

    Output:
        rho_dot_kms_vals (np.ndarray): 
            Array of ring intercept radial velocities in km/s.
        phi_rl_dot_kms_vals (np.ndarray): 
            Array of ring intercept azimuthal velocties in km/s.

    Note:
        [1] np.gradient assumes constant time spacing. From np.gradient
            documentation: The gradient is computed using second order 
            accurate central differences in the interior points and either 
            first or second order accurate one-sides (forward or backwards) 
            differences at the boundaries. The returned gradient hence has 
            the same shape as the input array.
            (https://docs.scipy.org/doc/numpy/reference/generated/numpy.gradient.html)
    """

    # Compute central differences using numpy gradient
    rho_dot_kms_vals = np.gradient(rho_km_vals, dt)
    phi_rl_rad_vals = np.radians(phi_rl_deg_vals)

    phi_rl_dot_kms_vals = rho_km_vals * np.gradient(phi_rl_rad_vals, dt)
    return rho_dot_kms_vals, phi_rl_dot_kms_vals
