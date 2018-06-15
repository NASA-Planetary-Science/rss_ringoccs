
'''

calc_occ_geometry.py

Purpose: Class Geom whose attributes are the occultation geometry parameters.

Revisions:
    2018 May 08 - jfong - copied from calc_geometry_v4.py
    2018 May 23 - jfong - make old methods into functions, make new
                          methods that calc the info in data_catalog
    2018 May 30 - jfong - add history attribute
    2018 Jun 01 - jfong - add calc_rho_corr_pole, calc_rho_corr_timing 
                          functions
                        - generalize planet and spacecraft
    2018 Jun 13 - jfong - use np.arctan2 instead of spice.vsep in 
                          xform_J2000_to_ring_plane_frame()
'''

import time
import os
import platform
import spiceypy as spice
import numpy as np
from ..tools.spm_to_et import spm_to_et
from ..tools.et_to_spm import et_to_spm
import pdb

#SPACECRAFT              = 'Cassini'
#PLANET                  = 'Saturn'
#SATURN_NAIF_ID          = 699
#SATURN_NAIF_RADII       = 'BODY699_RADII'

class Geometry(object):

    """This is an object that calculates occultation geometry from event time.
    """


    def __init__(self, rsr_inst, planet, spacecraft, kernels, pt_per_sec=1.):
        """This calculates occultation geometry as attributes.

        Args:
            rsr_inst (obj): Header of RSR file containing SPM and F_sky_hz
            kernels (list(str)): List of SPICE kernels that cover event time.
            planet 
            spacecraft
            pt_per_sec (int, optional): Number of points calculated per second for all attributes.
                This will default to 1.
        """

        self.kernels = kernels
        self.input_vars = {
                "rsr_inst": None, #rsr_inst.history,
                "kernels": kernels
                }
        self.input_kwds = {
                "pt_per_sec": pt_per_sec
                }
        self.history = self.write_hist_dict()
        self.rsr_inst = rsr_inst


        # Extract information from rsr instance
        year    = rsr_inst.year
        doy        = rsr_inst.doy
        dsn        = rsr_inst.dsn
        band    = rsr_inst.band
        sample_rate_khz = rsr_inst.sample_rate_khz
        rsr_file = rsr_inst.rsr_file
        spm_raw = rsr_inst.spm_vals
        (f_spm, f_sky) = rsr_inst.get_f_sky_pred()
        
        # Create new spm array with defined points per second        
        spm_start = spm_raw[0]
        spm_end = spm_raw[-1]
        step = 1./pt_per_sec
        t_oet_spm_vals = np.arange(spm_start, spm_end, step)
        t_oet_et_vals = spm_to_et(t_oet_spm_vals, doy, year, kernels=kernels)

        # Interpolate to get sky frequency
        F_sky_hz_vals = np.interp(t_oet_spm_vals, f_spm, f_sky)

        # Load kernels
        load_kernels(kernels)

        # Calculate Saturn center to ring intercept vector
        rho_vec_vals, t_ret_et_vals = calc_rho_vec_km(t_oet_et_vals, planet,
                spacecraft, dsn, kernels=kernels)

        
        rho_km_vals = []
        for rho_vec in rho_vec_vals:
            rho_km_vals.append(spice.vnorm(rho_vec))

        # Calculate spacecraft event time
        t_set_et_vals = calc_set_et(t_oet_et_vals, spacecraft, dsn)


        # Retrieve Saturn pole unit vector
        nhat_p = get_pole(t_set_et_vals[0], planet)

        phi_rl_deg_vals, phi_ora_deg_vals = calc_phi_deg(t_oet_et_vals,
                rho_vec_vals, spacecraft, dsn, nhat_p)

        D_km_vals = calc_D_km(t_ret_et_vals, t_set_et_vals)

        B_deg_vals = calc_B_deg(t_oet_et_vals, spacecraft, dsn, nhat_p)

        F_km_vals = calc_F_km(D_km_vals, F_sky_hz_vals, B_deg_vals, 
                phi_ora_deg_vals)

        t_ret_spm_vals = et_to_spm(t_ret_et_vals)
        t_set_spm_vals = et_to_spm(t_set_et_vals)

        rho_dot_kms_vals, phi_rl_dot_kms_vals = calc_rip_velocity(rho_km_vals,
                phi_rl_deg_vals, step)

        R_sc_km_vals, R_sc_dot_kms_vals = calc_sc_state(t_set_et_vals, 
                spacecraft, planet, dsn, nhat_p)

        R_imp_km_vals = calc_impact_radius_km(R_sc_km_vals, t_set_et_vals,
                spacecraft, dsn, nhat_p)
        alt_deg_vals = calc_altitude_deg(t_oet_et_vals, spacecraft, dsn)


        self.t_oet_spm_vals = np.asarray(t_oet_spm_vals)
        self.t_ret_spm_vals = np.asarray(t_ret_spm_vals)
        self.t_set_spm_vals = np.asarray(t_set_spm_vals)
        self.rho_km_vals = np.asarray(rho_km_vals)
        self.phi_rl_deg_vals = np.asarray(phi_rl_deg_vals)
        self.phi_ora_deg_vals = np.asarray(phi_ora_deg_vals)
        self.D_km_vals = np.asarray(D_km_vals)
        self.B_deg_vals = np.asarray(B_deg_vals)
        self.rho_dot_kms_vals = np.asarray(rho_dot_kms_vals)
        self.phi_rl_dot_kms_vals = np.asarray(phi_rl_dot_kms_vals)
        self.F_km_vals = np.asarray(F_km_vals)
        self.R_imp_km_vals = np.asarray(R_imp_km_vals)
        self.rx_km_vals = np.stack(R_sc_km_vals)[:,0]
        self.ry_km_vals = np.stack(R_sc_km_vals)[:,1]
        self.rz_km_vals = np.stack(R_sc_km_vals)[:,2]
        self.vx_kms_vals = np.stack(R_sc_dot_kms_vals)[:,0]
        self.vy_kms_vals = np.stack(R_sc_dot_kms_vals)[:,1]
        self.vz_kms_vals = np.stack(R_sc_dot_kms_vals)[:,2]

#        self.nhat_p = nhat_p
#        self.t_oet_et_vals = t_oet_et_vals
#        self.t_set_et_vals = t_set_et_vals
#        self.t_ret_et_vals = t_ret_et_vals
        self.alt_deg_vals = alt_deg_vals
#
#        self.get_saturn_occ_times()
#        a, b, c = self.get_dsn_info()
    def get_planet_naif_id(self, planet):
        return spice.bodn2c(planet)


    def write_hist_dict(self):
        """This creates a history dictionary.

        Returns:
            geo_hist (dict): Dictionary with "user name", "host name",
                    "run date", "python version", "operating system",
                    "source file", "input variables", and "input keywords".
        """
        user_name = os.getlogin() 
        host_name = os.uname()[1]
        run_date = time.ctime() + ' ' + time.tzname[0] 
        python_version = platform.python_version() 
        operating_system = os.uname()[0]          
        src_file = __file__ 

        geo_hist = {
                           "user name": user_name,
                           "host name": host_name,
                            "run date": run_date,
                      "python version": python_version,
                    "operating system": operating_system,
                         "source file": src_file,
                     "input variables": self.input_vars,
                      "input keywords": self.input_kwds                                                }
        return geo_hist

    def get_rev(self):
        """This gets the rev number of the occultation event, col 2 of
            data catalog.

        Returns:
            rev (int): Rev number of occultation event.
        """
        # TODO

        return None

    def get_year(self):
        """This gets the year of the occultation event, col 3 of data catalog.

        Returns:
            year (int): Year of occultation event.
        """
        return self.rsr_inst.year

    def get_doy(self):
        """This gets the day of year of occultation event, col 4 of data
            catalog.

        Returns:
            doy (int): Day of year of occultation event.
        """
        return self.rsr_inst.doy

    def get_utc_start(self):
        doy = self.rsr_inst.doy
        year = self.rsr_inst.year
        spm_start = self.t_oet_spm_vals[0]
        et_start = spm_to_et(spm_start, doy, year, kernels=kernels)

        format_str = 'C'
        prec = 3
        utc_start = spice.et2uc(et_start, format_str, prec)
        return utc_start

    def get_utc_end(self):
        doy = self.rsr_inst.doy
        year = self.rsr_inst.year
        spm_end = self.t_oet_spm_vals[-1]
        et_end = spm_to_et(spm_start, doy, year, kernels=kernels)

        format_str = 'C'
        prec = 3
        utc_end = spice.et2uc(et_start, format_str, prec)
        return utc_end

    def get_dsn_info(self):
        idnum = self.rsr_inst.dsn[-2:]
        id_list = ['34', '35', '36', '43', '45',
                '14','15','24','25','26',
                '54', '55', '63', '65',
                '84']
        cmplx_list = ['G', 'G', 'G', 'G', 'G',
                    'M', 'M', 'M', 'M',
                    'C', 'C', 'C', 'C', 'C',
                    'ML']
        aperture_list = ['70','34', '34', '34', '34',
                    '34', '34', '70', '34',
                    '34', '34', '34', '70', '34',
                    '35']
        ind = id_list.index(idnum)
        cmplx = cmplx_list[ind]
        aperture = aperture_list[ind]

    
        return idnum, cmplx, aperture

    def get_freq_band(self):
        return self.rsr_inst.band

    def get_rsr_file(self):
        return self.rsr_inst.rsr_file


    def get_sampling_rate(self):
        return self.sampling_rate_khz

    def get_spm_start(self):
        return self.t_oet_spm_vals[0]

    def get_spm_end(self):
        return self.t_oet_spm_vals[-1]

    def get_rho_start(self):
        return self.rho_km_vals[0]

    def get_rho_end(self):
        return self.rho_km_vals[-1]


    def get_min_alt(self):
        return min(self.altitude_deg)

    def get_max_alt_info(self):
        maxalt = max(self.altitude_deg_vals)
        ind = np.where(self.altitude_deg_vals == maxalt)
        spm = self.t_oet_spm_vals[ind]
        rho = self.rho_km_vals[ind]
        return maxalt, spm, rho

    def get_B_start(self):
        et_start = self.t_oet_et_vals[0]
        dsn = self.dsn
        nhat_p = self.nhat_p

        B_deg_start = calc_B_deg(et_start, dsn, nhat_p)

        return B_deg_start
    def get_saturn_occ_times(self):
        et_vals = self.t_set_et_vals
        dsn = self.dsn

        npts = len(et_vals)
        
        abcorr = 'CN'

        # Get original body radii
        dim, BODY699_RADII = spice.bodvrd('SATURN', 'RADII', 3)
        
        # Add atmosphere + ionosphere to radii
        #    use only ionosphere radii for a conservative estimate
#        BODY699_ATM = [rho699 + 500.0 for rho699 in BODY699_RADII]
        BODY699_ION = [rho699 + 5000.0 for rho699 in BODY699_RADII]

        # Use atmospheric radii in kernel pool
        spice.pdpool('BODY699_RADII', BODY699_ION)

        et_blocked = []
        for n in range(npts):
            et = et_vals[n]
            occ_code_ion = spice.occult('Saturn', 'ELLIPSOID', 
                    'IAU_SATURN', 'Cassini', 'POINT', ' ', abcorr, dsn, et)
            if occ_code_ion != 0:
                et_blocked.append(et)

        spm_blocked = et_to_spm(et_blocked)


        return spm_blocked
    
def calc_rho_km(et_vals, planet, spacecraft, dsn, kernels=None):
    rho_vec, t_ret = calc_rho_vec_km(et_vals, planet, spacecraft, dsn,
            kernels=kernels)
    rho_km = []
    for vec in rho_vec:
        rho_km.append(spice.vnorm(vec))

    return rho_km
    
def calc_rho_corr_ncfiv_eq4(rho_km, rho_dot_kms, delta_t, alpha): #, rmin, rmax):
    r0 = 100000.
#    term1 = rho + rho_dot*delt
    term1 = rho_dot_kms * delta_t
    term2 = alpha * ((rho_km - r0)/1000.)
    drho = term1 - term2
    return drho

def calc_rho_corr_timing(rho_dot_kms, delta_t):
    return rho_dot_kms * delta_t
    
def calc_rho_corr_pole(rho_km_ori, et, dsn, planet, spacecraft,  new_kernels):
    rho_vec_km_new, t_ret_et_new = calc_rho_vec_km(et, planet, spacecraft, dsn, kernels=new_kernels)
    rho_km_new = []
    for rho_vec in rho_vec_km_new:
        rho_km_new.append(spice.vnorm(rho_vec))

    rho_corr_km = rho_km_new - rho_km_ori
    return rho_corr_km

def load_kernels(kernels):
    """This clears the kernel pool and loads a new set of kernels.

    Args:
        kernels (list(str)): List of paths to NAIF kernels.
    
    Returns:
        None
    """
    # Clear kernels
    spice.kclear()

    # Load kernels
    for kernel in kernels:
        spice.furnsh(kernel)

    return None

def calc_altitude_deg(et_vals, spacecraft, dsn, kernels=None):
    
    npts = len(et_vals)
    elev_deg_vals = np.zeros(npts)
    
    if kernels:
        load_kernels(kernels)

    for n in range(npts):
        et = et_vals[n]
        # Compute dsn to spacecraft position vector in J2000
        ref = 'J2000'
        abcorr = 'CN'
        ptarg1, ltime1 = spice.spkpos(spacecraft, et, ref, abcorr, dsn)

        # Compute Earth to dsn position vector in J2000
        abcorr = 'NONE'
        planet = 'EARTH'

        ptarg2, ltime2 = spice.spkpos(dsn, et, ref, abcorr, planet)

        # Calculate elevation as the complement to the angle between
        #    ptarg1 (dsn->spacecraft) and ptarg2 (earth->dsn)
        elev_deg_vals[n] = 90. - spice.vsep(ptarg1, ptarg2)*spice.dpr()




    return elev_deg_vals

def get_pole(et, planet, kernels=None):
    """This calculates unit vector in pole direction from kernel constants.

    Args:
        et (float64): Ephemeris seconds past J2000

    Returns:
        nhat_p (np.ndarray): 1x3 unit vector in pole direction 
    """
    if kernels:
        load_kernels(kernels)

    # Retrieve right ascension and declination from kernel pool
    #     Note: the quadratic terms for pole direction are typically zero,
    #    but are retained here for consistency with PCK file format
    #    definitions.
    bodynm = planet
    item = 'POLE_RA'
    maxn = 3
    dim1, pole_RA = spice.bodvrd(bodynm, item, maxn)

    bodynm = planet
    item = 'POLE_DEC'
    maxn = 3
    dim2, pole_DEC = spice.bodvrd(bodynm, item, maxn)

    dt_centuries = et / (spice.spd() * 365.25*100.)
    RAP = (pole_RA[0] + dt_centuries*pole_RA[1]
            + dt_centuries**2*pole_RA[2])
    DEP = (pole_DEC[0] + dt_centuries*pole_DEC[1]
            + dt_centuries**2*pole_DEC[2])

    # Convert to rectangular coordinates
    inrange = 1.
    re = RAP * spice.rpd()
    dec = DEP * spice.rpd()
    nhat_p = spice.radrec(inrange, re, dec)

    return nhat_p

def calc_set_et(t_oet_et_vals, spacecraft, dsn, kernels=None):
    """This calculates the spacecraft event time.

    Args:
        t_oet_et_vals (np.ndarray): Array of Earth-received times in ET sec.
        dsn(str): ID of Deep Space Network receiver
    
    Returns:
        t_set_et_vals (np.ndarray): Array of spacecraft event times in ET sec.
    """

    npts = len(t_oet_et_vals)
    t_set_et_vals = np.zeros(npts)
    
    for n in range(npts):
        t_oet_et = t_oet_et_vals[n]
    
        # Compute ltime, the one-way light travel time between Cassini
        #    and dsn station
        targ    = spacecraft
        ref = 'J2000'
        abcorr  = 'CN'
        et = t_oet_et
        obs = dsn
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)
    
        nhat_sc2dsn = spice.vhat(starg)
    
        t_set_et_vals[n] = t_oet_et - ltime

    return t_set_et_vals

def calc_rho_vec_km(t_oet_et_vals, planet, spacecraft, dsn, kernels=None):
    """This calculates rho_vec, the position vector (in J2000 frame) from
            Saturn center to the intersection of a light ray from 
            Cassini to the input DSN and the Saturn ring plane.

    Args:
        t_oet_et_vals(np.ndarray): Array of earth-received times in ET sec.
        dsn(str): ID of Deep Space Network receiver

    Returns:
        rho_vec_vals(list(np.ndarray)): list of 3xN arrays of the 
            position vector from Saturn center to the intersection 
            of a light ray from Cassini to the input DSN and the
            Saturn ring plane.
        t_ret_et_vals(np.ndarray): Ring event time in ET seconds.
    """
    if kernels:
        load_kernels(kernels)

    planet_id = spice.bodn2c(planet)

    npts = len(t_oet_et_vals)
    rho_vec_vals = [] #np.zeros((3, npts))
    t_ret_et_vals = np.zeros(npts)

    # Replace Saturn radii in kernel pool to values that
    #    represent a flat ellisoid, which will represent the
    #     Saturn ring plane.
    body = planet_id
    item = 'RADII'
    dim = 3
    radii = spice.bodvcd(body, item, dim)
    
    new_radii = [1.e6, 1.e6, 1.e-5]

    # #SATURN_NAIF_RADII       = 'BODY699_RADII'
    planet_naif_radii = 'BODY'+str(planet_id)+'_RADII'

    name = planet_naif_radii
    dvals = new_radii
    spice.pdpool(name, dvals)
        

    for n in range(npts):
        t_oet_et = t_oet_et_vals[n]
        # Compute starg, the spacecraft position relative to DSN 
        #     at station at Earth-received time to Cassini, with light 
        #    correction, and light travel time 
        targ    = spacecraft
        ref = 'J2000'
        abcorr  = 'CN'
        et = t_oet_et
        obs = dsn
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        nhat_sc2dsn = spice.vhat(starg)

        # Compute the intersection of vector with the ring plane,
        #     the time epoch in ET seconds of intersection (aka ring event
        #     time), and vector from observer to intercept point.
        iau_planet = 'IAU_'+ planet.upper()

        method   = 'Ellipsoid'
        target   = planet
        et          = t_oet_et
#        fixref   = 'IAU_SATURN'
        fixref   = iau_planet
        abcorr   = 'CN'
        obsrvr   = dsn
        dref     = 'J2000'
        dvec     = nhat_sc2dsn
        spoint, trgepc, srfvec = spice.sincpt(method, target,
                t_oet_et, fixref, abcorr, obsrvr, dref, dvec)
        
        t_ret_et_vals[n] = trgepc
    
        # Convert ring plane intercept to J2000 frame
        frame_from = fixref
        frame_to = dref
        etfrom = trgepc
        etto = t_oet_et
        xform = spice.pxfrm2(frame_from, frame_to, etfrom, etto)
        
        rho_vec = spice.mxv(xform, spoint)
        rho_vec_vals.append(rho_vec)
        
    # Restore old valuoes of RADII to kernel pool
#    name = SATURN_NAIF_RADII
    name = planet_naif_radii
    dvals = radii[1].tolist()
    spice.pdpool(name, dvals)
        
    return rho_vec_vals, t_ret_et_vals

def calc_phi_deg(t_oet_et_vals, rho_km_vec_vals, spacecraft, dsn, nhat_p, kernels=None):
    """This calculates observed ring azimuth and ring longitude.

    Args:
        t_oet_et_vals (np.ndarray): Array of Earth-received time in ET sec.
        rho_vec_vals (np.ndarray): Nx3 array of ring intercept point position
            vector in km.

    Returns:
        phi_rl_deg_vals (np.ndarray): Array of ring longitude in degrees.
        phi_ora_deg_vals (np.ndarray): Array of observed ring azimuth 
            in degrees.
    """
    if kernels:
        load_kernels(kernels)

    npts = len(t_oet_et_vals)
    phi_rl_deg_vals = np.zeros(npts)
    phi_ora_deg_vals = np.zeros(npts)


    for n in range(npts):
        t_oet_et = t_oet_et_vals[n]

        rho_km_vec = rho_km_vec_vals[n]

        # Compute dsn to Cassini position vector with light correction
        targ    = dsn
        et = t_oet_et
        ref = 'J2000'
        abcorr  = 'CN'
        obs = spacecraft
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        # Rotate vector so that xy plane is in the ring plane
        zaxis = [0., 0., 1.]
        axdef = nhat_p
        indexa = 3
        plndef = zaxis
        indexp = 2
        mm13 = spice.twovec(axdef, indexa, plndef, indexp)
        vec_RL = spice.mxv(mm13, rho_km_vec)
        
        radius_vec_RL, RA_vec_RL, DEC_vec_RL = spice.recrad(vec_RL)

        phi_rl_deg = RA_vec_RL * spice.dpr()
        phi_rl_deg_vals[n] = phi_rl_deg

        vec_ORA = spice.mxv(mm13, starg)
        vec_ORA[2] = 0.
        radius_vec_ORA, RA_vec_ORA, DEC_vec_ORA = spice.recrad(vec_ORA)

        # Note: phi_ora_deg differs from the MTR86 definition by 180 deg
        phi_ora_deg = (phi_rl_deg - RA_vec_ORA*spice.dpr() + 720.) % 360.
        phi_ora_deg_vals[n] = phi_ora_deg




    return phi_rl_deg_vals, phi_ora_deg_vals

def calc_D_km(t1, t2):
    """This calculates the light distance between two input times.

    Args:
        t1(np.ndarray): Array of time in seconds. 
        t2(np.ndarray): Array of time in seconds.

    Returns:
        D_km_vals(np.ndarray): Array of light distance in kilometers.
    """
    npts = len(t1)
    D_km_vals = np.zeros(npts)
    for n in range(npts):
        dt = abs(t1[n] - t2[n])
        D_km = dt * spice.clight()
        D_km_vals[n] = D_km

    return D_km_vals

def calc_B_deg(et_vals, spacecraft, dsn, nhat_p, kernels=None):
    """This calculates ring opening angle.

    Args:
        t_oet_vals (np.ndarray): Array of ET seconds
        dsn(str): Deep space station receiver ID

    Returns:
        B_deg_vals (np.ndarray): Array of ring opening angle in degrees.
    """
    if kernels:
        load_kernels(kernels)

    npts = len(et_vals)
    B_deg_vals = np.zeros(npts)


    for n in range(npts):
        # Compute Cassini to dsn position vector
        targ = dsn
        et = et_vals[n]
        ref = 'J2000'
        abcorr = 'CN'
        obs = spacecraft
        starg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)

        # Calculate B as the complement to the angle
        #     made by the Saturn pole vector and the 
        #    Cassini to dsn vector
        v1 = starg
        v2 = nhat_p
        B_rad = (np.pi/2.) - spice.vsep(v1, v2)

        B_deg_vals[n] = B_rad * spice.dpr()

    return B_deg_vals

def calc_F_km(D_km_vals, F_sky_hz_vals, B_deg_vals, phi_ora_deg_vals):
    """This calculates the Fresnel scale (see MTR1986 Eq 6).

    Args:
        D_km_vals (np.ndarray): Array of light distance in kilometers.
        F_sky_hz_vals (np.ndarray): Array of downlink sinusoidal signal
            frequency at the front-end of dsn station.
        B_deg_vals (np.ndarray): Array of ring opening angle in degrees.
        phi_ora_deg_vals (np.ndarray): Array of observed ring azimuth 
            in degrees.
    
    Returns:
        F_km_vals (np.ndarray): Array of Fresnel scale in km.
    """
    npts = len(D_km_vals)
    F_km_vals = np.zeros(npts)

    for n in range(npts):
        lambda_sky = spice.clight() / F_sky_hz_vals[n]
        phi_ORA_rad = spice.rpd() * phi_ora_deg_vals[n]
        B_rad = spice.rpd() * B_deg_vals[n]
        D_km = D_km_vals[n]

        F_km = np.sqrt((0.5* lambda_sky * D_km * 
                        (1-(np.cos(B_rad))**2* 
                        (np.sin(phi_ORA_rad))**2))/ 
                        (np.sin(B_rad))**2)
        F_km_vals[n] = F_km
    return F_km_vals

def xform_J2000_to_ring_plane_frame(vec, et, spacecraft, dsn, nhat_p, kernels=None):
    """This transforms a vector from J2000 frame to Saturn ring plane frame.

    Args:
        vec(np.ndarray): 3-element vector in J2000 frame
        et(float64): Ephemeris time corresponding to given vec
        dsn(str): ID of Deep Space Station receiver
    
    Returns:
        out_vec(np.ndarray): 3-element vector in Saturn ring plane frame
    """
    if kernels:
        load_kernels(kernels)
    # Define planetocentric frame as z-axis in Saturn pole direction
    # and x-axis in Cassini to dsn direction
    zaxis = [0., 0., 1.]
    xaxis = [1., 0., 0.]

    # Compute Cassini (at et) to dsn state vector (at et+lighttime) 
    targ = dsn
    et = et
    ref = 'J2000' 
    abcorr = 'XCN'
    obs = spacecraft
    starg1, ltime1 = spice.spkezr(targ, et, ref, abcorr, obs)

    # Rotate z-axis to Saturn pole direction
    axdef = nhat_p
    indexa = 3
    plndef = zaxis
    index = 2
    rotmat_z = spice.twovec(axdef, indexa, plndef, index)
    R_sat2sc_km_z = spice.mxv(rotmat_z, vec)
    R_sc2dsn_km_z = spice.mxv(rotmat_z, starg1[0:3])

    # Rotate x-axis to Cassini to dsn direction
    nhat_sc2dsn = R_sc2dsn_km_z/np.linalg.norm(R_sc2dsn_km_z)
    nhat_sc2dsn[2] = 0.

#    rot_angle = -spice.vsep(nhat_sc2dsn, xaxis)
    rot_angle = np.arctan2(nhat_sc2dsn[1], nhat_sc2dsn[0])
    rot_mat_sc2dsn = spice.rotate(rot_angle, 3)

    out_vec = spice.mxv(rot_mat_sc2dsn, R_sat2sc_km_z)
    return out_vec


def calc_sc_state(t_set_et_vals, spacecraft, planet, dsn, nhat_p, kernels=None):
    """This calculates Cassini state vector in a planetocentric frame.

    Args:
        t_set_et_vals(np.ndarray): Array of spacecraft event time in ET sec
        dsn(str): ID of Deep Space Station receiver
    
    Returns:
        R_sc_km_vals(list(np.ndarray)): List of Cassini position vectors in
            km.
        R_sc_dot_kms_vals(list(np.ndarray)): List of Cassini velocity 
            vectors in km/s.
    """
    if kernels:
        load_kernels(kernels)

    npts = len(t_set_et_vals)
    R_sc_km_vals = []
    R_sc_dot_kms_vals = []
    for n in range(npts):
        t_set_et = t_set_et_vals[n]
        # Saturn to Cassini state vector, with no light-time correction
        targ = spacecraft
        et = t_set_et
        ref = 'J2000'
        abcorr = 'NONE'
        obs = planet
        starg0, ltime0 = spice.spkezr(targ, et, ref, abcorr, obs)
        
        R_sat2sc_km = starg0[0:3]
        R_sat2sc_dot_kms = starg0[3:6]

        # Transform vectors to Saturn ring plane frame
        R_sc_km = xform_J2000_to_ring_plane_frame(
                R_sat2sc_km, t_set_et, spacecraft, dsn, nhat_p)
        R_sc_dot_kms = xform_J2000_to_ring_plane_frame(
                R_sat2sc_dot_kms, t_set_et, spacecraft, dsn, nhat_p)

        R_sc_km_vals.append(R_sc_km)
        R_sc_dot_kms_vals.append(R_sc_dot_kms)

    return R_sc_km_vals, R_sc_dot_kms_vals
def calc_rip_velocity(rho_km_vals, phi_rl_deg_vals, dt):
    """This calculates the ring intercept point radial and azimuthal
            velocity.

    Args:
        rho_km_vals(np.ndarray): Array of ring intercept points in km.
        phi_rl_deg_vals(np.ndarray): Array of ring longitudes in deg.
        dt(float64): Time spacing in sec between points, must be constant.
    
    Returns:
        rho_dot_kms_vals(np.ndarray): Array of ring intercept radial 
            velocities in km/s.
        phi_rl_dot_kms_vals(np.ndarray): Array of ring intercept 
            azimuthal velocities in km/s.
    """

    rho_dot_kms_vals = np.gradient(rho_km_vals, dt)
    phi_rl_rad_vals = phi_rl_deg_vals * spice.rpd()

    phi_rl_dot_kms_vals = rho_km_vals * np.gradient(phi_rl_rad_vals, dt)
    return rho_dot_kms_vals, phi_rl_dot_kms_vals

def calc_impact_radius_km(R_sc_km_vals, et_vals, spacecraft, dsn, nhat_p, kernels=None):
    """This calculates the impact radius of in km.

    Args:
        et_vals(np.ndarray): Array of Earth-received time in ET sec.
        R_sc_km_vals(list(np.ndarray)): List of 3x1 arrays of spacecraft
            position vector in planetocentric frame at given et_vals.
        dsn(str): ID of Deep Space Station receiver
    
    Returns:
        R_imp_km_vals(np.ndarray): Array of impact radius in km.
    """

    if kernels:
        load_kernels(kernels)

    npts = len(et_vals)
    R_imp_km_vals = np.zeros(npts)
    for n in range(npts):
        R_sc_km = R_sc_km_vals[n]
        # Compute Cassini to dsn position vector with no light correction
        targ = dsn
        et = et_vals[n]
        ref = 'J2000' 
        abcorr = 'XCN'
        obs = spacecraft
        starg1, ltime1 = spice.spkezr(targ, et, ref, abcorr, obs)

        # Transform vector to Saturn ring plane frame
        R_sc2dsn_km_pcf = xform_J2000_to_ring_plane_frame(
                starg1[0:3], et, spacecraft, dsn, nhat_p)

        lindir = R_sc2dsn_km_pcf
        linpt = R_sc_km
        point = [0., 0., 0.]
        pnear, distance = spice.nplnpt(linpt, lindir, point)

        R_imp_km_vals[n] = distance

    return R_imp_km_vals
