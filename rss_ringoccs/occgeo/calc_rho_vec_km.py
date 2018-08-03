"""
calc_rho_vec_km.py

Purpose: Calculate the position vector of the intersection of the spacecraft-
         -to-observer ray with the planet ring plane, from the planet center.

Revisions:
    2018 Jul 09 - jfong - original
"""

import numpy as np
import spiceypy as spice


def calc_rho_vec_km(et_vals, planet, spacecraft, dsn, kernels=None,
        verbose=False):
    """
    This calculates the position vector of the ring intercept point from the
    planet center.

    Args:
        et_vals (np.ndarray): Array of earth-received times in ET sec
        planet (str): Name of planet -- must be compatible with NAIF.
        spacecraft (str): Name of spacecraft -- must be compatible with NAIF
        dsn (str): DSN observing station ID -- must be compatible with NAIF
        kernels (list): List of NAIF kernels, including path.
        verbose (bool): Default False prints nothing, while setting as
                              True prints steps

    Output:
        rho_vec_km_vals (list): List of 3xN np.ndarrays of the planet center to
                                ring intercept point position vector
        t_ret_et_vals (np.ndarray): Array of ring event times in ET seconds.

    References:
        [1] Ring intercept point calculation using a dynamical frame:
            page 19 of https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/27_derived_quant.pdf
    """

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernels)

    planet_id = spice.bodn2c(planet)

    npts = len(et_vals)
    rho_vec_km_vals = []
    t_ret_et_vals = np.zeros(npts)

    # Replace Saturn radii in kernel pool with values that represent a
    #   flat ellipsoid, which we will use as the ring plane
    body = planet_id
    item = 'RADII'
    dim = 3
    radii = spice.bodvcd(body, item, dim)

    new_radii = [1.e6, 1.e6, 1.e-5]

    # Construct naif radii code in form: 'BODY699_RADII'
    planet_naif_radii = 'BODY'+str(planet_id)+'_RADII'

    name = planet_naif_radii
    dvals = new_radii
    spice.pdpool(name, dvals)

    for n in range(npts):
        et = et_vals[n]

        # Compute spacecraft position relative to dsn
        targ = spacecraft
        ref = 'J2000'
        abcorr = 'CN'
        obs = dsn
        starg, ltime, = spice.spkpos(targ, et, ref, abcorr, obs)

        nhat_sc2dsn = spice.vhat(starg)

        # Compute intersection of vector with ring plane and time epoch
        #   of intersection in ET secs (ring event time)
        iau_planet = 'IAU_'+planet.upper()

        method = 'Ellipsoid'
        target = planet
        fixref = iau_planet
        abcorr = 'CN'
        obsrvr = dsn
        dref = 'J2000'
        dvec = nhat_sc2dsn
        spoint, trgepc, srfvec = spice.sincpt(method, target, et,
                fixref, abcorr, obsrvr, dref, dvec)

        t_ret_et_vals[n] = trgepc

        # Convert ring plane intercept to J2000 frame
        frame_from = fixref
        frame_to = dref
        etfrom = trgepc
        etto = et
        xform = spice.pxfrm2(frame_from, frame_to, etfrom, etto)

        rho_vec_km = spice.mxv(xform, spoint)
        rho_vec_km_vals.append(rho_vec_km)

    # Restore old original values of RADII to kernel pool
    name = planet_naif_radii
    dvals = radii[1].tolist()
    spice.pdpool(name, dvals)

    return rho_vec_km_vals, t_ret_et_vals
