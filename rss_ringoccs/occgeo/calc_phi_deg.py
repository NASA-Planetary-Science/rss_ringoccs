"""

calc_phi_deg.py

Purpose: Calculate observed ring azimuth, the angle measured from the
         direction of a photon heading to the observer to the direction
         of a local radial vector (projected into ring plane and measured in
         the prograde direction). This function also calculates the
         intertial longitude relative to the prime meridian.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_phi_deg(et_vals, rho_vec_km_vals, spacecraft, dsn, nhat_p,
        kernels=None):
    """
    This calculates observed ring azimuth and ring longitude.

    Args:
        et_vals (np.ndarray): 
            Array of earth-received time in ET seconds.
        rho_vec_km_vals (np.ndarray): 
            Nx3 array of ring intercept position vectors in km.
        spacecraft (str): 
            Name of spacecraft -- must be compatible with NAIF.
        dsn (str): 
            DSN observing station ID -- must be compatible with NAIF
        nhat_p (np.ndarray): 
            1x3 array unit vector in planet pole direction.
        kernels (list): 
            List of NAIF kernels, including path

    Outputs:
        phi_rl_deg_vals (np.ndarray): 
            Array of inertial longitude in degrees.

        phi_ora_deg_vals (np.ndarray): 
            Array of observed ring azimuth in degrees.

    Note:
        [1] phi_ora_deg differns from the MTR1986 definition by 180 degrees.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    npts = len(et_vals)
    phi_rl_deg_vals = np.zeros(npts)
    phi_ora_deg_vals = np.zeros(npts)

    for n in range(npts):
        et = et_vals[n]

        rho_vec_km = rho_vec_km_vals[n]

        # Compute DSN to spacecraft position vector with light correction
        targ = dsn
        ref = 'J2000'
        abcorr = 'CN'
        obs = spacecraft
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        # Rotate vector so that xy plane is in the ring plane
        zaxis = [0., 0., 1.]
        axdef = nhat_p
        indexa = 3
        plndef = zaxis
        indexp = 2

        rotmat = spice.twovec(axdef, indexa, plndef, indexp)
        vec_rl = spice.mxv(rotmat, rho_vec_km)

        # Convert coordinates to RA and DEC for inertial longitude
        radius_vec_rl, ra_vec_rl, dec_vec_rl = spice.recrad(vec_rl)

        phi_rl_deg = ra_vec_rl * spice.dpr()
        phi_rl_deg_vals[n] = phi_rl_deg

        # Calculate observed ring azimuth by rotating pm to direction of
        #   photon heading to observer
        vec_ora = spice.mxv(rotmat, starg)

        # Project vector onto ring plane
        vec_ora[2] = 0

        # Convert coordinates to RA and DEC for ORA
        radius_vec_ora, ra_vec_ora, dec_vec_ora = spice.recrad(vec_ora)

        phi_ora_deg_vals[n] = (phi_rl_deg - ra_vec_ora*spice.dpr()
                + 720.) % 360.

    return phi_rl_deg_vals, phi_ora_deg_vals
