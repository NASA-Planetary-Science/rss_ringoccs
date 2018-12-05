"""

rss_occgeo.py

Purpose: Calculate occultation geometry for RSS ring events.

Revisions:
    2018 Jul 10 - jfong - original
    2018 Jul 11 - jfong - remove write_hist method and use
                          write_history_dict() from tools/
    2018 Jul 12 - jfong - add planetary_occultation_flag and
                          ring_occultation_direction and
                          ring_profile_direction and
                          revolution_number attributes and methods
    2018 Jul 16 - jfong - add beta_vals attribute
    2018 Jul 26 - jfong - remove write_history_dict() import
                        - add __write_history_dict method
    2018 Jul 27 - jfong - remove __write_history_dict method and use 
                          write_history_dict() again, which now requires
                          source file as an additional input
    2018 Aug 02 - jfong - remove methods pertaining to data catalog entries
                        - remove planet_occ_flag, rev_num, ring_occ_dir,
                          rsr_inst attributes
                        - remove rsr_file and sample_rate_khz attribute 
                          extraction from rsr_inst
    2018 Sep 20 - jfong - add write_file kwarg
                        - set rev_info from rsr_inst as an attribute
                        - add prof_dir to rev_info dict attribute
    2018 Oct 15 - jfong - add freespace_km attribute using find_gaps()
    2018 Oct 17 - jfong - add ionos_occ_spm_vals, atmos_occ_spm_vals attr
    2018 Oct 23 - jfong - update freespace_km attr to freespace_spm
                        - add ingress_spm, egress_spm attributes
    2018 Nov 05 - jfong - use SET for atmosphere calculation
    2018 Nov 19 - jfong - add chord verification
                        - replace splrep, splev for interp1d for rho interp

"""
from ..tools.spm_to_et import spm_to_et
from ..tools.et_to_spm import et_to_spm
from ..tools.write_output_files import write_output_files
from ..tools.write_history_dict import write_history_dict
#from ..tools.write_output_files import write_output_files
#from ..tools.write_history_dict_v2 import write_history_dict_v2

#from ..tools.pds3_geo_series import write_geo_series
from .calc_elevation_deg import calc_elevation_deg
from .get_pole import get_pole
from .calc_rho_vec_km import calc_rho_vec_km
from .calc_set_et import calc_set_et
from .calc_phi_deg import calc_phi_deg
from .calc_D_km import calc_D_km
from .calc_B_deg import calc_B_deg
from .calc_F_km import calc_F_km
from .xform_j2k_to_pcf import xform_j2k_to_pcf
from .calc_sc_state import calc_sc_state
from .calc_rip_velocity import calc_rip_velocity
from .calc_impact_radius_km import calc_impact_radius_km
from .get_planet_occ_times import get_planet_occ_times
from .calc_beta import calc_beta
from .calc_beta import calc_B_eff_deg
from .get_freespace import get_freespace
from scipy.interpolate import splrep
from scipy.interpolate import splev


import spiceypy as spice
import numpy as np
import sys
import pdb
import inspect

class Geometry(object):

    """
    This is an object that calculates occultation geometry needed for
    diffraction reconstruction as well as other relevant geometry parameters.
    """
    def __init__(self, rsr_inst, planet, spacecraft, kernels, pt_per_sec=1.,
            verbose=False, write_file=True):
        """This calculates occultation geometry as attributes.

        Args:
            rsr_inst (class): 
                Instance of RSRReader class.
            kernels (list): 
                List of NAIF kernels, including path.
            planet (str): 
                Planet name -- must be compatible with NAIF.
            spacecraft (str):
                Spacecraft name -- must be compatible with NAIF
            pt_per_sec (float64, optional): 
                Number of points calculated per second for all geometry
                calculations. This will default to 1.
            verbose (bool):
                Boolean for whether processing steps are printed.

        Attributes:
            t_oet_spm_vals (np.ndarray):
                Observed event time in seconds past midnight.
            t_ret_spm_vals (np.ndarray):
                Ring event time in seconds past midnight.
            t_set_spm_vals (np.ndarray):    
                Spacecraft event time in seconds past midnight.
            rho_km_vals (np.ndarray):       
                Distance in km from the center of Saturn to ring intercept
                point.
            phi_rl_deg_vals (np.ndarray):
                Ring longitude (inertial longitude) in degrees.
            phi_ora_deg_vals (np.ndarray):
                Observed ring azimuth in degrees.
            D_km_vals (np.ndarray):         
                Spacecraft to ring intercept point distance in km.
            B_deg_vals (np.ndarray):        
                Ring opening angle in degrees.
            rho_dot_kms_vals (np.ndarray): 
                Ring intercept radial velocity in km/s.
            phi_rl_dot_kms_vals (np.ndarray):
                Ring intercept azimuthal velocity in km/s.
            F_km_vals (np.ndarray): 
                Fresnel scale in km.
            R_imp_km_vals (np.ndarray): 
                Impact radius in km.
            rx_km_vals (np.ndarray):
                x-component of spacecraft position in a planetocentric frame,
                in km.
            ry_km_vals (np.ndarray):        
                y-component of spacecraft position in a planetocentric frame,
                in km.
            rz_km_vals (np.ndarray):        
                z-component of spacecraft position in a planetocentric frame,
                in km.
            vx_kms_vals (np.ndarray):       
                x-component of spacecraft velocity in a planetocentric frame,
                in km/s.
            vy_kms_vals (np.ndarray):       
                y-component of spacecraft velocity in a planetocentric frame,
                in km/s
            vz_kms_vals (np.ndarray):       
                z-component of spacecraft velocity in a planetocentric frame,
                in km/s
            elev_deg_vals (np.ndarray):
                Elevation angle in degrees.
            history (dict):                 
                Dictionary of processing history.
            naif_toolkit_version (str):     
                NAIF toolkit version used (e.g., "V.N0066").
            B_eff_deg_vals (np.ndarray):    
                Effective ring opening angle in deg.
            beta_vals (np.ndarray):         
                Optical depth enhancement factor.
            kernels (list):
                List of NAIF kernels, including path.
            frequency_band (str):   
                Letter initial of frequency band of the radio signal sent by 
                the spacecraft (e.g.,  'X', 'K', 'S')

        Dependencies:
            [1] numpy
            [2] rss_ringoccs.occgeo
            [3] rss_ringoccs.tools
            [4] spiceypy
            [5] sys

        Notes:
            [1] kernels list must include:
                1) spacecraft ephemeris kernel
                2) planetary constants kernel
                3) leapseconds kernel
                4) planet and lunar ephemeris kernel
                5) earth stations kernel
                6) earth rotation and constants kernel

        Warnings:
            [1] This code has only been tested for planet='Saturn'.
        """
        #pdb.set_trace()
        #ex=inspect.getargvalues(inspect.currentframe())[3]

        #pdb.set_trace()
        # write input history

        if verbose:
            print('Calculating occultation geometry...')

        if verbose:
            print('\tChecking for valid inputs...')

        if type(planet) != str:
            sys.exit('WARNING (Geometry): Input planet is NOT a string!')

        if type(spacecraft) != str:
            sys.exit('WARNING (Geometry): Input spacecraft is NOT a string!')

        if not isinstance(pt_per_sec, (int, float)):
            sys.exit('WARNING (Geometry): Input pt_per_sec is NOT an int or '
                    + 'float!')

        if verbose:
            print('\tExtracting information from rsr file...')

        # Extract information from rsr instance
        try:
            year = rsr_inst.year
            doy = rsr_inst.doy
            dsn = rsr_inst.dsn
            band = rsr_inst.band
            self.__spm_full = rsr_inst.spm_vals
            #self.__et_vals = spm_to_et(self.__spm_full, doy, year, kernels=kernels)
            spm_start = rsr_inst.spm_vals[0]
            spm_end = rsr_inst.spm_vals[-1]
            rsr_hist = rsr_inst.history
            (f_spm, f_sky) = rsr_inst.get_f_sky_pred()

            if verbose:
                print('\tSetting rev info attribute...')

            # Setting rev info attribute
            self.rev_info = rsr_inst.rev_info


        except AttributeError:
            sys.exit('WARNING (Geometry): Input RSR instance does '
                    + 'not have valid attributes!')

        # Create new spm array with defined points per second
        step = 1./pt_per_sec
        t_oet_spm_vals = np.arange(spm_start, spm_end, step)
        t_oet_et_vals = spm_to_et(t_oet_spm_vals, doy, year, kernels=kernels)
        #self.__et_vals = np.arange(t_oet_et_vals[0], t_oet_et_vals[-1],
        #        0.1)

        # Interpolate to get sky frequency
        f_sky_hz_vals = np.interp(t_oet_spm_vals, f_spm, f_sky)

        # Calculate Saturn center to ring intercept vector
        if verbose:
            print('\tCalculating ring intercept event time and vector...')
        rho_vec_vals, t_ret_et_vals = calc_rho_vec_km(t_oet_et_vals, planet,
                spacecraft, dsn, kernels=kernels)

        rho_km_vals = [spice.vnorm(vec) for vec in rho_vec_vals]

        # Calculate spacecraft event time
        if verbose:
            print('\tCalculating spacecraft event time...')
        t_set_et_vals = calc_set_et(t_oet_et_vals, spacecraft, dsn)

        # Retrieve Saturn pole unit vector
        if verbose:
            print('\tRetrieving Saturn pole unit vector...')
        nhat_p = get_pole(t_set_et_vals[0], planet)

        if verbose:
            print('\tCalculating ring longitude and ring azimuth...')

        # Calculated ring longitude and observed ring azimuth
        phi_rl_deg_vals, phi_ora_deg_vals = calc_phi_deg(t_oet_et_vals,
                rho_vec_vals, spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating distance from spacecraft to ring intercept '
                + 'point...')
        # Calculate distance from spacecraft to ring intercept point
        D_km_vals = calc_D_km(t_ret_et_vals, t_set_et_vals)

        if verbose:
            print('\tCalculating ring opening angle...')

        # Calculate ring opening angle
        B_deg_vals = calc_B_deg(t_oet_et_vals, spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating Fresnel scale...')

        # Calculate Fresnel scale
        F_km_vals = calc_F_km(D_km_vals, f_sky_hz_vals, B_deg_vals,
                phi_ora_deg_vals)

        t_ret_spm_vals = et_to_spm(t_ret_et_vals)
        t_set_spm_vals = et_to_spm(t_set_et_vals)

        if verbose:
            print('\tCalculating ring intercept velocities...')
        # Calculate ring intercept velocities
        rho_dot_kms_vals, phi_rl_dot_kms_vals = calc_rip_velocity(rho_km_vals,
                phi_rl_deg_vals, step)

        if verbose:
            print('\tCalculating spacecraft state vector...')
        # Calculate spacecraft state vector
        R_sc_km_vals, R_sc_dot_kms_vals = calc_sc_state(t_set_et_vals,
                spacecraft, planet, dsn, nhat_p)

        if verbose:
            print('\tCalculating impact radius...')
        # Calculate impact radius
        R_imp_km_vals = calc_impact_radius_km(R_sc_km_vals, t_set_et_vals,
                spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating elevation angle...')

        # Calculate target angle above the horizon
        elev_deg_vals = calc_elevation_deg(t_oet_et_vals, spacecraft, dsn)

        # Calculate beta
        if verbose:
            print('\tCalculating beta...')
        beta_vals = calc_beta(B_deg_vals, phi_ora_deg_vals)
        B_eff_deg_vals = calc_B_eff_deg(B_deg_vals, phi_ora_deg_vals)

        # Set attributes, first block contains the inputs to *GEO.TAB
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
        self.rx_km_vals = np.stack(R_sc_km_vals)[:, 0]
        self.ry_km_vals = np.stack(R_sc_km_vals)[:, 1]
        self.rz_km_vals = np.stack(R_sc_km_vals)[:, 2]
        self.vx_kms_vals = np.stack(R_sc_dot_kms_vals)[:, 0]
        self.vy_kms_vals = np.stack(R_sc_dot_kms_vals)[:, 1]
        self.vz_kms_vals = np.stack(R_sc_dot_kms_vals)[:, 2]

        self.kernels = kernels
        self.frequency_band = band
        self.elev_deg_vals = np.asarray(elev_deg_vals)
        self.naif_toolkit_version = self.__get_naif_version()
        self.beta_vals = np.asarray(beta_vals)
        self.B_eff_deg_vals = np.asarray(B_eff_deg_vals)



        # Calculate when signal passes atmosphere + ionosphere
        ionos_occ_et_vals = get_planet_occ_times(t_oet_et_vals, dsn,
                planet, spacecraft, height_above=5000.)
        self.ionos_occ_spm_vals = et_to_spm(ionos_occ_et_vals)

        atmos_occ_et_vals = get_planet_occ_times(t_oet_et_vals, dsn,
                planet, spacecraft, height_above=500.)
        self.atmos_occ_spm_vals = et_to_spm(atmos_occ_et_vals)
        # Write processing history dictionary attribute
        if verbose:
            print('\tWriting history dictionary...')

        input_vars = {
                "rsr_inst": rsr_hist,
                "kernels": kernels,
                "planet": planet,
                "spacecraft": spacecraft
                }
        input_kwds = {
                "pt_per_sec": pt_per_sec
                }
        self.history = write_history_dict(input_vars, input_kwds, __file__)

        # Add prof_dir entry to rev_info dict
        self.rev_info['prof_dir'] = self.get_profile_dir()

        if self.rev_info['prof_dir'] == '"BOTH"':
            self.split_ind = self.get_chord_ind()
        else:
            self.split_ind = None

        self.freespace_km, self.freespace_spm = get_freespace(
                t_ret_spm_vals, year, doy, rho_km_vals,
                phi_rl_deg_vals, t_oet_spm_vals, self.atmos_occ_spm_vals,
                split_ind=self.split_ind, kernels=kernels)

        # Write output data and label file if set
        if write_file:
            write_output_files(self)

    def __get_naif_version(self):
        """
        This returns the NAIF toolkit version used.
        """
        naif_ver = spice.tkvrsn('TOOLKIT')
        return ('"V.' + naif_ver.split('_')[-1] + '"')

    def get_profile_dir(self):
        """
        This returns the profile direction observed.
        """
        rho = self.rho_km_vals
        dr_start = rho[1] - rho[0]
        dr_end = rho[-1] - rho[-2]

        if (dr_start < 0) and (dr_end < 0):
            prof_dir = '"INGRESS"'
        elif (dr_start > 0) and (dr_end > 0):
            prof_dir = '"EGRESS"'
        else:
            # check if actually a chord occultation
            prof_dir = self.verify_chord()
        return prof_dir

    def verify_chord(self):
        split_ind = self.get_chord_ind()

        t_oet_ing = self.t_oet_spm_vals[:split_ind]
        t_oet_egr = self.t_oet_spm_vals[split_ind:]

        n_ing = len(t_oet_ing)
        n_egr = len(t_oet_egr)


        atmos_blocked_spm = list(self.atmos_occ_spm_vals)
        
        ing_blocked = [x for x in t_oet_ing if x in atmos_blocked_spm]
        egr_blocked = [x for x in t_oet_egr if x in atmos_blocked_spm]

        if len(ing_blocked) == n_ing:
            prof_dir = '"EGRESS"'
        elif len(egr_blocked) == n_egr:
            prof_dir = '"INGRESS"'
        else:
            prof_dir = '"BOTH"'
    
        return prof_dir

    def get_chord_ind(self):


        # add 2 for the index of the first element after the sign change
        ind = np.argwhere(np.diff(np.sign(self.rho_dot_kms_vals)))

        if len(ind) > 1:
            print('WARNING! ring radius changes direction twice!')
            pdb.set_trace()
        elif len(ind) == 0:
            ind = None
        else:
            ind = int(ind[0]) + 2

        return ind

