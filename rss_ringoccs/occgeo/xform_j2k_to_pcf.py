"""

xform_j2k_to_pcf.py

Purpose: Transform a vector in J2000 frame to ring plane frame (or
         planetocentric frame). The xaxis is defined to be parallel to the
         projection of the line-of-sight from the spacecraft at SET to Earth
         at OET on the ring plane and the zaxis is defined by the direction
         of the planet's pole. The yaxis is defined by the right-hand rule.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def xform_j2k_to_pcf(vec, et, spacecraft, dsn, nhat_p, kernels=None):
    """
    This transforms a vector in J2000 frame to planet ring plane frame.

    Args:
        vec (np.ndarray):
            3-element vector in J2000 frame
        et (float64): 
            ET in seconds corresponding to input vec
        dsn (str): 
            DSN observing station ID -- must be compatible with NAIF.
        nhat_p (np.ndarray): 
            1x3 array unit vector in planet pole direction.
        kernels (list): 
            List of NAIF kernels, including path.

    Output:
        out_vec (np.ndarray): 
            3-element vector in planet ring plane frame.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Compute Cassini (at et) to dsn position vector (at et+ltime)
    targ = dsn
    ref = 'J2000'
    abcorr = 'XCN'
    obs = spacecraft
    starg1, ltime1 = spice.spkpos(targ, et, ref, abcorr, obs)

    # Rotate z-axis to planet pole direction
    zaxis = [0., 0., 1.]

    axdef = nhat_p
    indexa = 3
    plndef = zaxis
    index = 2
    rotmat_z = spice.twovec(axdef, indexa, plndef, index)
    vec_z = spice.mxv(rotmat_z, vec)
    R_sc2dsn_km_z = spice.mxv(rotmat_z, starg1)

    # Rotate xaxis to Cassini to dsn direction
    nhat_sc2dsn = R_sc2dsn_km_z/np.linalg.norm(R_sc2dsn_km_z)
    nhat_sc2dsn[2] = 0.

    rot_angle = np.arctan2(nhat_sc2dsn[1], nhat_sc2dsn[0])
    rot_mat_sc2dsn = spice.rotate(rot_angle, 3)

    out_vec = spice.mxv(rot_mat_sc2dsn, vec_z)

    return out_vec
