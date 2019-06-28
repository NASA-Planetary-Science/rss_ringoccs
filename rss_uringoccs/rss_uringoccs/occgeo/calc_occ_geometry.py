"""
:Purpose:

    Functions for calculating occultation geometry parameters listed in
    *GEO.LBL file from Archived_Cassini_RSS_RingOccs_2018/ and other
    useful geometry information, such as free-space regions and planetary
    occultation times.

:Dependencies:
    #. spiceypy
    #. numpy

"""
import spiceypy as spice
import numpy as np
from scipy.interpolate import interp1d

def calc_B_deg(et_vals, spacecraft, dsn, nhat_p, kernels=None,
        ref='J2000'):
    """
    This calculates ring opening angle, or the observed ring elevation,
    as the complement to the angle made by the planet pole vector and
    the spacecraft to DSN vector

    Arguments:
        :et_vals (*np.ndarray*): Array of earth-received times in
            ephemeris seconds.
        :dsn (*str*): DSN observing station ID -- must be compatible with NAIF.
        :nhat_p (*np.ndarray*): 1x3 array unit vector in planet pole direction.

    Keyword Arguments:
        :kernels (*str* or *list*): List of NAIF kernels, including path.

    Returns:
        :B_deg_vals (*np.ndarray*): Array of ring opening angle in degrees.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    if isinstance(et_vals, float):
        npts = 1
        et_vals = [et_vals]
    else:
        npts = len(et_vals)

    if npts == 0:
        return []

    B_deg_vals = np.zeros(npts)

    for n in range(npts):
        # Compute spacecraft to dsn position vector
        targ = dsn
        et = et_vals[n]
        #ref = 'J2000'
        abcorr = 'CN'
        obs = spacecraft
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        # Calculate B as the complement to the angle made by the
        #   planet pole vector and the spacecraft to DSN vector
        v1 = starg
        v2 = nhat_p
        B_rad = (np.pi/2.) - spice.vsep(v1, v2)

        B_deg_vals[n] = B_rad * spice.dpr()

    return B_deg_vals

def calc_D_km(t1, t2):
    """
    This calculates the light distance between two input times, in km.

    Args:
        :t1 (*np.ndarray*): Array of time in seconds.
        :t2 (*np.ndarray*): Array of time in seconds.

    Returns:
        :D_km_vals (*np.ndarray*): Array of light distance in km.
    """
    D_km_vals = abs(t1-t2) * spice.clight()
    return D_km_vals

def calc_F_km(D_km_vals, f_sky_hz_vals, B_deg_vals, phi_ora_deg_vals):
    """
    This calculates the Fresnel scale using Eq. 6 of [MTR1986]_.

    Arguments:
        :D_km_vals (*np.ndarray*): Array of spacecraft to ring intercept
            point distances in km.
        :f_sky_hz_vals (*np.ndarray*): Array of downlink sinusoidal signal
            frequency at the front-end  of observing dsn station, in Hz.
        :B_deg_vals (*np.ndarray*): Array of ring opening angle in degrees.
        :phi_ora_deg_vals (*np.ndarray*): Array of observed ring azimuth
            in degrees.

    Returns:
        :F_km_vals (*np.ndarray*): Array of Fresnel scale in km.

    Notes:
        #. diffcorr uses an independently-calculated Fresnel scale
        #. Reference: [MTR1986]_ Equation 6

    """

    lambda_sky = spice.clight() / f_sky_hz_vals
    phi_ora_rad = spice.rpd() * phi_ora_deg_vals
    B_rad = spice.rpd() * B_deg_vals
    D_km = D_km_vals

    # Eq. 6 of MTR1986
    F_km_vals = np.sqrt((0.5 * lambda_sky * D_km * (1 - (np.cos(B_rad))**2
        * (np.sin(phi_ora_rad))**2)) / (np.sin(B_rad))**2)

    return F_km_vals

def calc_B_eff_deg(B_deg, phi_ora_deg):
    """
    This calculates the effective ring opening angle in degrees.

    Arguments:
        :B_deg (*float* or *np.ndarray*): Ring opening angle in degrees.
        :phi_ora_deg (*float* or *np.ndarray*): Observed ring azimuth in
            degrees.
    Returns:
        :B_eff_deg (*float* or *np.ndarray*): Effective ring opening angle
            in degrees.

    Notes:
        #. Reference: [GRESH86]_ Eq. 16
    """
    B_rad = np.radians(B_deg)
    phi_ora_rad = np.radians(phi_ora_deg)
    B_eff_deg = np.arctan2(np.tan(B_rad), np.cos(phi_ora_rad))
    return np.asarray(B_eff_deg)

def calc_beta(B_deg, phi_ora_deg):
    """
    This calculates the optical depth enhancement factor.

    Arguments:
        :B_deg (*float* or *np.ndarray*): Ring opening angle in degrees.
        :phi_ora_deg (*float* or *np.ndarray*): Observed ring azimuth in
            degrees.

    Returns:
        :beta (*np.ndarray*): Optical depth enhancement factor.
    """
    B_eff_deg = calc_B_eff_deg(B_deg, phi_ora_deg)
    beta = 1./(np.tan(B_eff_deg))
    return np.asarray(beta)


def calc_elevation_deg(et_vals, target, obs, kernels=None):
    """
    Calculate the elevation of a target above the horizon for a given observer.

    Arguments:
        :et_vals (*np.ndarray*): Array of observed event times in ET sec.
        :target (*str*): Target name -- must be compatible with NAIF. This
            will typically be spacecraft or planet name.
        :obs (*str*): Observer name -- must be compatible with NAIF. This
            will typically be observing dsn station name. Observer is assumed
            to be on Earth.
        :kernels (*str* or *list*): List of NAIF kernels, including path.

    Returns:
        :elev_deg_vals (*np.ndarray*): Array of elevation angles in degrees.
    """
    npts = len(et_vals)
    elev_deg_vals = np.zeros(npts)

    # Load kernels
    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    for n in range(npts):
        et = et_vals[n]

        # Compute observer to target position vector in J2000
        #   with light-correction
        ref = obs + '_TOPO'
        abcorr = 'CN'
        ptarg1, ltime1 = spice.spkpos(target, et, ref, abcorr, obs)

        # Compute Earth to observer position vector in J2000
        #   without light correction
        abcorr = 'NONE'
        planet = 'EARTH'

        ptarg2, ltime2 = spice.spkpos(obs, et, ref, abcorr, planet)

        # Calculate elevation as the complement to the angle
        #   between ptarg1 (obs->target) and ptarg2 (Earth->obs)
        elev_deg_vals[n] = 90. - spice.vsep(ptarg1, ptarg2)*spice.dpr()

    return elev_deg_vals

def calc_impact_radius_km(R_sc_km_vals, et_vals, spacecraft, dsn, nhat_p,
        ref='J2000', kernels=None):
    """
    This calculates the closest approach of the spacecraft signal to the
    planet defined as a sphere.

    Arguments:
        :R_sc_km_vals (*list*): List of 3-element arrays of spacecraft
            position vector in planetocentric frame at input et_vals.
        :et_vals (*np.ndarray*): Array of Earth-received times in ephemeris
            seconds.
        :spacecraft (*str*): Spacecraft name
        :dsn (*str*): DSN observing station ID
        :nhat_p (*np.ndarray*): 1x3 array unit vector in planet pole direction.

    Keyword Arguments:
        :kernels (*list*): List of NAIF kernels, including path.

    Returns:
        :R_imp_km_vals (*np.ndarray*): Array of impact radius in km.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    npts = len(et_vals)
    R_imp_km_vals = np.zeros(npts)

    for n in range(npts):
        R_sc_km = R_sc_km_vals[n]
        et = et_vals[n]

        # Compute spacecraft to dsn position vector in J2000 frame,
        #   at et+ltime
        targ = dsn
        #ref = 'J2000'
        abcorr = 'XCN'
        obs = spacecraft
        starg1, ltime1 = spice.spkpos(targ, et, ref, abcorr, obs)

        # Transform vector to planetocentric frame
        R_sc2dsn_km_pcf = xform_j2k_to_pcf(starg1, et, spacecraft,
                dsn, nhat_p)

        # Calculate distance from saturn center to point of closest approach
        lindir = R_sc2dsn_km_pcf
        linpt = R_sc_km
        point = [0., 0., 0.]
        pnear, distance = spice.nplnpt(linpt, lindir, point)

        R_imp_km_vals[n] = distance
    return R_imp_km_vals

def calc_phi_deg(et_vals, rho_vec_km_vals, spacecraft, dsn, nhat_p, ref='J2000',
        kernels=None):
    """
    This calculates observed ring azimuth and ring longitude.

    Arguments:
        :et_vals (*np.ndarray*): Array of earth-received time in ET seconds.
        :rho_vec_km_vals (*np.ndarray*): Nx3 array of ring intercept position
            vectors in km.
        :spacecraft (*str*): Name of spacecraft
        :dsn (*str*): DSN observing station ID
        :nhat_p (*np.ndarray*): 1x3 array unit vector in planet pole direction.

    Keyword Arguments:
        :kernels (*str* or *list*): List of NAIF kernels, including path

    Returns:
        :phi_rl_deg_vals (*np.ndarray*): Array of inertial longitude in degrees.
        :phi_ora_deg_vals (*np.ndarray*): Array of observed ring azimuth
            in degrees.

    Notes:
        #. phi_ora_deg differs from the [MTR1986]_ definition by 180 degrees.
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
        #ref = 'J2000'
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

def calc_rho_km(et_vals, planet, spacecraft, dsn, kernels=None,
        ring_frame=None):
    """
    Calculate the distance between Saturn center to ring intercept point.

    Arguments:
        :et_vals (*np.ndarray*): Array of observed event times in ET sec.
        :planet (*str*): Planet name
        :spacecraft (*str*): Spacecraft name
        :dsn (*str*): DSN observing station ID

    Keyword Arguments:
        :kernels (*str* or *list*): List of NAIF kernels, including path.

    Returns:
        :rho_km_vals (*np.ndarray*): Array of ring intercept points in km.
    """

    # Calculate rho vector
    rho_vec_km, t_ret_et = calc_rho_vec_km(et_vals, planet, spacecraft, dsn,
            kernels=kernels, ring_frame=ring_frame)

    # Compute magnitude of each vector
    rho_km_vals = [spice.vnorm(vec) for vec in rho_vec_km]

    return np.asarray(rho_km_vals)

def calc_rho_vec_km(et_vals, planet, spacecraft, dsn, ref='J2000', kernels=None,
        verbose=False, ring_frame=None):
    """
    This calculates the position vector of the ring intercept point from the
    planet center in J2000 frame.

    Arguments:
        :et_vals (*np.ndarray*): Array of earth-received times in ET sec
        :planet (*str*): Name of planet
        :spacecraft (*str*): Name of spacecraft
        :dsn (*str*): DSN observing station ID

    Keyword Arguments:
        :kernels (*str* or *list*): Path to NAIF kernels
        :verbose (*bool*): Option for printing processing steps

    Output:
        :rho_vec_km_vals (*list*): List of 3xN np.ndarrays of the planet
            center to ring intercept point position vector in J2000 frame
        :t_ret_et_vals (*np.ndarray*): Array of ring event times in ET seconds.

    References:
        #. Ring intercept point calculation using a dynamical frame.
            See [NAIF]_ page 19.

    """

    if kernels:
        spice.kclear()
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
    if ring_frame is None:
        iau_planet = 'IAU_'+planet.upper()
        frame = iau_planet
    else:
        frame = ring_frame
    #print(frame)


    for n in range(npts):
        et = et_vals[n]

        # Compute spacecraft position relative to dsn
        targ = spacecraft
        #ref = 'J2000'
        abcorr = 'CN'
        obs = dsn
        starg, ltime, = spice.spkpos(targ, et, ref, abcorr, obs)

        nhat_sc2dsn = spice.vhat(starg)

        # Compute intersection of vector with ring plane and time epoch
        #   of intersection in ET secs (ring event time)

        method = 'Ellipsoid'
        target = planet
        #fixref = iau_planet
        fixref = frame
        abcorr = 'CN'
        obsrvr = dsn
        #dref = 'J2000'
        dref = ref
        dvec = nhat_sc2dsn
        spoint, trgepc, srfvec = spice.sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)

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

def calc_rip_velocity(rho_km_vals, phi_rl_deg_vals, dt):
    """
    This calculates the ring intercept point radial and azimuthal velocity.

    Arguments:
        :rho_km_vals (*np.ndarray*): Array of ring intercept points in km.
        :phi_rl_deg_vals (*np.ndarray*): Array of ring longitudes in degrees.
        :dt (*float*): Constant time spacing between points.

    Returns:
        :rho_dot_kms_vals (*np.ndarray*): Array of ring intercept radial
            velocities in km/s.
        :phi_rl_dot_kms_vals (*np.ndarray*): Array of ring intercept azimuthal
            velocties in km/s.

    Note:
        #. np.gradient assumes constant time spacing. From np.gradient
            documentation: The gradient is computed using second order
            accurate central differences in the interior points and either
            first or second order accurate one-sides (forward or backwards)
            differences at the boundaries. The returned gradient hence has
            the same shape as the input array.
            .. (https://docs.scipy.org/doc/numpy/reference/generated/numpy.gradient.html)
    """

    # Compute central differences using numpy gradient
    rho_dot_kms_vals = np.gradient(rho_km_vals, dt)
    phi_rl_rad_vals = np.radians(phi_rl_deg_vals)

    phi_rl_dot_kms_vals = rho_km_vals * np.gradient(phi_rl_rad_vals, dt)
    return rho_dot_kms_vals, phi_rl_dot_kms_vals

def calc_sc_state(et_vals, spacecraft, planet, dsn, nhat_p, ref='J2000',
        kernels=None):
    """
    This calculates spacecraft state vector in a planetocentric frame.

    Arguments:
        :et_vals (*np.ndarray*): Array of spacecraft event times in ET seconds.
        :spacecraft (*str*): Spacecraft name
        :planet (*str*): Planet name
        :dsn (*str*): Deep Space Network observing station ID
        :nhat_p (*np.ndarray*): 1x3 array unit vector in planet pole direction.


    Keyword Arguments:
        :kernels (*str* or *list*): Path to NAIF kernel(s)

    Returns:
        :R_sc_km_vals (*list*): List of Nx3 np.ndarrays of spacecraft position
            vector in km in planetocentric frame
        :R_sc_dot_kms_vals (*list*): List of Nx3 np.ndarrays of spacecraft
            velocity vector in km/s.

    Notes:
        #.  Saturn planetocentric frame is defined by x-axis parallel to
            projection of spacecraft-to-Earth line-of-sight, z-axis in
            direction of Saturn's pole.
    """

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernel)

    npts = len(et_vals)
    R_sc_km_vals = []
    R_sc_dot_kms_vals = []
    for n in range(npts):
        et = et_vals[n]

        # Compute planet to spacecraft state vector in J2000 frame,
        #   with no light-time correction
        targ = spacecraft
        #ref = 'J2000'
        abcorr = 'NONE'
        obs = planet
        starg0, ltime0 = spice.spkezr(targ, et, ref, abcorr, obs)

        R_sat2sc_km = starg0[0:3]
        R_sat2sc_dot_kms = starg0[3:6]

        # Transform vectors to planetocentric frame
        R_sc_km = xform_j2k_to_pcf(R_sat2sc_km, et, spacecraft, dsn, nhat_p)
        R_sc_dot_kms = xform_j2k_to_pcf(R_sat2sc_dot_kms, et,
                spacecraft, dsn, nhat_p)

        R_sc_km_vals.append(R_sc_km)
        R_sc_dot_kms_vals.append(R_sc_dot_kms)

    return R_sc_km_vals, R_sc_dot_kms_vals

def calc_set_et(et_vals, spacecraft, dsn, kernels=None):
    """
    This calculates the time at which photons left the spacecraft, in ET sec.

    Arguments:
        :et_vals (*np.ndarray*): Array of earth-received times in ET seconds.
        :spacecraft (*str*):
        :dsn (*str*): Deep Space Network observing station ID

    Keyword Arguments:
        :kernels (*str* or *list): Path to NAIF kernels

    Returns:
        :t_set_et_vals (*np.ndarray*): Array of spacecraft event times
            in ET sec.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Convert spacecraft and dsn strings to NAIF integer codes
    sc_code = spice.bodn2c(spacecraft)
    dsn_code = spice.bodn2c(dsn)

    # Use light travel time to calculate spacecraft event time
    t_set_et_vals = []
    for et in et_vals:
        set_et, ltime = spice.ltime(et, dsn_code, "<-", sc_code)
        t_set_et_vals.append(set_et)

    return np.asarray(t_set_et_vals)

def get_start_jd(year, doy):
    """
    Purpose:
        Get the start of a day in Julian date times.

    Arguments:
        :year (*str*): Year
        :doy (*str*): Day of year

    Returns:
        :start_jd (*float*): Julian date time of the start of a day
    """

    et = spice.utc2et(str(year)+'-'+(str(doy)).zfill(3)+'T00:00:00')
    jd = spice.et2utc(et, 'J', 10)
    start_jd = float(jd[3:])
    return start_jd

def split_chord_arr(t_ret_spm_vals, t_oet_spm_vals,
        atmos_occ_spm_vals, phi_rl_deg_vals, rho_km_vals, ind, profdir):
    """
    Purpose:
        Return array of only ingress or egress portion of a chord occultation.

    Arguments:
        :t_ret_spm_vals (*np.ndarray*): Ring event times in SPM
        :t_oet_spm_vals (*np.ndarray*): Observed event times in SPM
        :atmos_occ_spm_vals (*np.ndarray*): SPM times of when spacecraft signal
            is blocked by planet atmosphere
        :phi_rl_deg_vals (*np.ndarray*): Inertial ring longitude in deg.
        :rho_km_vals (*np.ndarray*): Ring intercept points in km
        :ind (*int*): Index of where ingress switches to egress
        :profdir (*str*): Profile direction to return ('"INGRESS"' or
            '"EGRESS"')

    Returns:
        :t_ret_spm_vals (*np.ndarray*): Ring event times in SPM of 'profdir'
            portion of occultation
        :t_oet_spm_vals (*np.ndarray*): Observed event times in SPM of 'profdir'
            portion of occultation
        :phi_rl_deg_vals (*np.ndarray*): Inertial ring longitude in deg of
            'profdir' portion of occultation
        :rho_km_vals (*np.ndarray*): Ring intercept points in km of 'profdir'
            portion of occultation
    """

    # Split egress portion
    if profdir == '"EGRESS"':
        t_ret_spm_split = t_ret_spm_vals[ind:]
        t_oet_spm_split = t_oet_spm_vals[ind:]
        phi_rl_deg_split = phi_rl_deg_vals[ind:]
        rho_km_split = rho_km_vals[ind:]

    # Split ingress portion
    if profdir == '"INGRESS"':
        t_ret_spm_split = t_ret_spm_vals[:ind]
        t_oet_spm_split = t_oet_spm_vals[:ind]
        phi_rl_deg_split = phi_rl_deg_vals[:ind]
        rho_km_split = rho_km_vals[:ind]

    # Remove any part of the occultation that is blocked by atmosphere
    t_ret_out, t_oet_out, phi_rl_out, rho_out = remove_blocked(
            t_oet_spm_split, atmos_occ_spm_vals, t_ret_spm_split,
            phi_rl_deg_split, rho_km_split)

    return t_ret_out, t_oet_out, phi_rl_out, rho_out

def remove_blocked(t_oet_spm_vals, atmos_occ_spm_vals, t_ret_spm_vals,
        phi_rl_deg_vals, rho_km_vals):
    """
    Purpose:
        Remove values that occur during times blocked by planet atmosphere.
    Arguments:
        :t_oet_spm_vals (*np.ndarray*): Observed event times in SPM
        :atmos_occ_spm_vals (*np.ndarray*): SPM times of when spacecraft signal
            is blocked by planet atmosphere
        :t_ret_spm_vals (*np.ndarray*): Ring event times in SPM
        :phi_rl_deg_vals (*np.ndarray*): Inertial ring longitude in deg.
        :rho_km_vals (*np.ndarray*): Ring intercept points in km

    Returns:
        :t_ret_spm_vals (*np.ndarray*): Ring event times in SPM, excluding
            atmospheric occultation times
        :t_oet_spm_vals (*np.ndarray*): Observed event times in SPM, excluding
            atmospheric occultation times
        :phi_rl_deg_vals (*np.ndarray*): Inertial ring longitude in deg,
            excluding atmospheric occultation times
        :rho_km_vals (*np.ndarray*): Ring intercept points in km, excluding
            atmospheric occultation times
    """
    t_oet_spm_vals = np.asarray(t_oet_spm_vals)
    atmos_occ_spm_vals = np.asarray(atmos_occ_spm_vals)
    t_ret_spm_vals = np.asarray(t_ret_spm_vals)
    phi_rl_deg_vals = np.asarray(phi_rl_deg_vals)
    rho_km_vals = np.asarray(rho_km_vals)

    npts = len(t_oet_spm_vals)
    mask_not_blocked = np.array([False for i in range(npts)])
    for i in range(npts):
        if t_oet_spm_vals[i] not in atmos_occ_spm_vals:
            mask_not_blocked[i] = True


    t_ret_out = t_ret_spm_vals[mask_not_blocked]
    t_oet_out = t_oet_spm_vals[mask_not_blocked]
    phi_rl_out = phi_rl_deg_vals[mask_not_blocked]
    rho_out = rho_km_vals[mask_not_blocked]

    return t_ret_out, t_oet_out, phi_rl_out, rho_out

def get_planet_occ_times(et_vals, obs, planet, spacecraft, height_above=500.,
        kernels=None):
    """
    Purpose:
        Return times when the spacecraft-to-observer ray is blocked by planet.

    Arguments:
        :et_vals (*np.ndarray*): Array of observed event times in ET sec.
        :obs (*str*): Observer name
        :planet (*str*): Planet name
        :spacecraft (*str*): Spacecraft name

    Keyword Arguments:
        :height_above (*float*): Height in km to be added to planet radius to
            account for the atmosphere
        :kernels (*str* or *list*): Path to NAIF kernels

    Returns:
        :et_blocked_vals (*np.ndarray*): Array of observed event times in ET

    Note:
        #. This was made to be generalizable to different planets, but has
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

def get_pole(et, planet, kernels=None):
    """
    Purpose:
        Calculate unit vector in pole direction from kernel constants.

    Arguments:
        :et (*float*): Ephemeris seconds past J2000
        :planet (*str*): Planet name

    Keyword Arguments:
        :kernels (*str* or *list*): Path to NAIF kernels

    Returns:
        :nhat_p (*np.ndarray*): 1x3 unit vector in pole direction.

    Note:
        #. Quadratic terms for pole direction are typically zero but
            are retained here for consistency with PCK file format definitions.
    """

    # Load kernels
    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Retrieve right ascension and declination from kernel pool
    bodynm = planet
    item = 'POLE_RA'
    maxn = 3
    dim1, pole_RA = spice.bodvrd(bodynm, item, maxn)

    bodynm = planet
    item = 'POLE_DEC'
    maxn = 3
    dim2, pole_DEC = spice.bodvrd(bodynm, item, maxn)

    # Calculate pole ra and dec using quadratic terms
    dt_centuries = et / (spice.spd()*365.25*100.)
    rap = (pole_RA[0] + dt_centuries*pole_RA[1]
            + dt_centuries**2*pole_RA[2])
    dep = (pole_DEC[0] + dt_centuries*pole_DEC[1]
            + dt_centuries**2*pole_DEC[2])

    # Convert to rectangular coordinates
    inrange = 1.
    re = rap * spice.rpd()
    dec = dep * spice.rpd()
    nhat_p = spice.radrec(inrange, re, dec)

    return nhat_p

def xform_j2k_to_pcf(vec, et, spacecraft, dsn, nhat_p, ref='J2000',
        kernels=None):
    """
    Purpose
        Transform vector in J2000 frame to planet ring plane frame.

    Arguments:
        :vec (*np.ndarray*): 3-element vector in J2000 frame
        :et (*float*): ET in seconds corresponding to input vec
        :dsn (*str*): DSN observing station ID
        :nhat_p (*np.ndarray*): 1x3 array unit vector in planet pole direction.

    Keyword Arguments:
        :kernels (*str* or *list*): Path to NAIF kernels

    Returns:
        :out_vec (*np.ndarray*): 3-element vector in planet ring plane frame.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # Compute Cassini (at et) to dsn position vector (at et+ltime)
    targ = dsn
    #ref = 'J2000'
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
