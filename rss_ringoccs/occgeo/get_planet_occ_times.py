"""

get_planet_occ_times.py

Purpose: Given earth-received time, get the subset of times at which the
         spacecraft-to-observer ray is blocked by the input planet.

Revisions:
    2018 Jul 09 - jfong - original
    2018 Jul 12 - jfong - generalize planet and spacecraft
    2018 Jul 31 - jfong - furnsh kernels if kwd provided
"""

import numpy as np
import spiceypy as spice


def get_planet_occ_times(et_vals, obs, planet, spacecraft, height_above=500.,
        kernels=None):
    """
    Get the times at which the spacecraft-to-observer ray is blocked by planet.

    Args:
        et_vals (np.ndarray):   
            Array of observed event times in ET sec.
        obs (str):              
            Observer name -- must be compatible with NAIF.
        planet (str):           
            Planet name -- must be compatible with NAIF
        spacecraft (str):       
            Spacecraft name -- must be compatible with NAIF
        height_above (float64): 
            Height in km to be added to planet radius to account for 
            the atmosphere. Default is 500km.
        kernels (list):         
            List of NAIF kernels, including path.

    Output:
        et_blocked_vals (np.ndarray): 
            Array of observed event times in ET sec at which the 
            spacecraft-to-obs ray is blocked by the input planet.

    Warnings:
        [1] This was made to be generalizable to different planets, but has
            been only tested with planet='Saturn'.

    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    npts = len(et_vals)
    et_blocked_vals = []

    planet_code = spice.bodn2c(planet)

    # Get planet radii from kernel pool
    dim, original_radii = spice.bodvrd(planet, 'RADII', 3)

    # Add height_above kwd to account for atmosphere
    new_radii = [radius + height_above for radius in original_radii]

    # Update kernel pool with new radii
    planet_naif_radii = 'BODY'+str(planet_code)+'_RADII'
    spice.pdpool(planet_naif_radii, new_radii)

    iau_planet = 'IAU_' + planet.upper()

    for n in range(npts):
        et = et_vals[n]
        # Determine occultation condition
        occ_code = spice.occult(planet, 'ELLIPSOID', iau_planet,
                spacecraft, 'POINT', ' ', 'CN', obs, et)
        # Append to et_blocked_vals if not partially or fully occulted
        if occ_code != 0:
            et_blocked_vals.append(et)

    return np.asarray(et_blocked_vals)
