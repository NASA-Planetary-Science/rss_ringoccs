"""
occgeo.py

:Purpose:

    Calculate occultation geometry for RSS ring events.

:Notes:
    #. kernels list must include:
        1) spacecraft ephemeris kernel
        2) planetary constants kernel
        3) leapseconds kernel
        4) planet and lunar ephemeris kernel
        5) earth stations kernel
        6) earth rotation and constants kernel
        7) topocentric frame kernel

:Dependencies:
    #. scipy.interpolate
    #. numpy
    #. spiceypy
    #. sys

"""
from ..tools.spm_to_et import spm_to_et
from ..tools.et_to_spm import et_to_spm
from ..tools.write_output_files import write_output_files
from ..tools.history import write_history_dict

import calc_occ_geometry as cog

from scipy.interpolate import splrep
from scipy.interpolate import splev

import spiceypy as spice
import numpy as np
import sys

class Geometry(object):

    """
    :Purpose:
        This is an object that calculates occultation geometry needed for
        diffraction reconstruction as well as other relevant geometry
        parameters.

    Arguments:
        :rsr_inst (*class*): Instance of RSRReader class.
        :kernels (*str* or *list*): List of NAIF kernels, including path.
        :planet (*str*): Planet name
        :spacecraft (*str*): Spacecraft name

    Keyword Arguments:
        :pt_per_sec (*float*): Number of points calculated per second
            for all geometry calculations.
        :verbose (*bool*): Boolean for whether processing steps are printed.
        :write_file (*bool*): Boolean for whether output *GEO.TAB and
            *GEO.LBL files will be created.

    Attributes:
        :t_oet_spm_vals (*np.ndarray*): Observed event time in seconds
            past midnight.
        :t_ret_spm_vals (*np.ndarray*): Ring event time in seconds past
            midnight.
        :t_set_spm_vals (*np.ndarray*): Spacecraft event time in seconds
            past midnight.
        :rho_km_vals (*np.ndarray*): Distance in km from the center of
            Saturn to ring intercept point.
        :phi_rl_deg_vals (*np.ndarray*): Ring longitude (inertial longitude)
            in degrees.
        :phi_ora_deg_vals (*np.ndarray*): Observed ring azimuth in degrees.
        :D_km_vals (*np.ndarray*): Spacecraft to ring intercept point
            distance in km.
        :B_deg_vals (*np.ndarray*): Ring opening angle in degrees.
        :rho_dot_kms_vals (*np.ndarray*): Ring intercept radial velocity
            in km/s.
        :phi_rl_dot_kms_vals (*np.ndarray*): Ring intercept azimuthal velocity
            in km/s.
        :F_km_vals (*np.ndarray*): Fresnel scale in km.
        :R_imp_km_vals (*np.ndarray*): Impact radius in km.
        :rx_km_vals (*np.ndarray*): x-component of spacecraft position in
            a planetocentric frame, in km.
        :ry_km_vals (*np.ndarray*): y-component of spacecraft position in
            a planetocentric frame, in km.
        :rz_km_vals (*np.ndarray*): z-component of spacecraft position in
            a planetocentric frame, in km.
        :vx_kms_vals (*np.ndarray*): x-component of spacecraft velocity in
            a planetocentric frame, in km/s.
        :vy_kms_vals (*np.ndarray*): y-component of spacecraft velocity in
            a planetocentric frame, in km/s
        :vz_kms_vals (*np.ndarray*): z-component of spacecraft velocity in
            a planetocentric frame, in km/s
        :elev_deg_vals (*np.ndarray*): Elevation angle in degrees.
        :kernels (*str* or *list*): List of NAIF kernels, including path.
        :rev_info (*dict*): RSR file specific info
        :history (*dict*): Dictionary of processing history.
        :naif_toolkit_version (*str*): NAIF toolkit version used
            (e.g., "V.N0066").
        :B_eff_deg_vals (*np.ndarray*): Effective ring opening angle in deg.
        :beta_vals (*np.ndarray*): Optical depth enhancement factor.
        :ionos_occ_spm_vals (*np.ndarray*): Array of seconds past midnight
            when the signal is occulted by the planet ionosphere, defined as
            5000km above an ellipsoid with radii from cpck file.
        :atmos_occ_spm_vals (*np.ndarray*): Array of seconds past midnight
            when the signal is occulted by the planet atmosphere, defined as
            500km above an ellipsoid with radii from cpck file.
        :freespace_spm (*list*): List of 2x1 lists of seconds past midnight
            values that define the inner and outer edge of a free-space
            gap in the ring system
        :freespace_km (*np.ndarray*): An array of 2x1 lists of km values
            that define the inner and outer edge of a free-space gap
            in the ring system


    """
    def __init__(self, rsr_inst, planet, spacecraft, kernels, pt_per_sec=1.,
            verbose=False, write_file=True):

        self.verbose = verbose

        if verbose:
            print('Calculating occultation geometry...')

        if verbose:
            print('\tChecking for valid inputs...')

        if type(planet) != str:
            raise ValueError('ERROR (Geometry): Input planet is NOT a string!')

        if type(spacecraft) != str:
            raise ValueError('ERROR (Geometry): Input spacecraft is NOT '
                            + 'a string!')

        if not isinstance(pt_per_sec, (int, float)):
            raise ValueError('ERROR (Geometry): Input pt_per_sec is NOT an int '
                                + 'or float!')

        if verbose:
            print('\tExtracting information from rsr file...')

        # Extract information from rsr instance
        year = rsr_inst.year
        doy = rsr_inst.doy
        dsn = rsr_inst.dsn
        band = rsr_inst.band
        spm_full = rsr_inst.spm_vals
        spm_start = rsr_inst.spm_vals[0]
        spm_end = rsr_inst.spm_vals[-1]
        rsr_hist = rsr_inst.history
        (f_spm, f_sky) = rsr_inst.get_f_sky_pred()

        self.rev_info = rsr_inst.rev_info

        # Create new spm array with defined points per second
        step = 1./pt_per_sec
        t_oet_spm_vals = np.arange(spm_start, spm_end, step)
        t_oet_et_vals = spm_to_et(t_oet_spm_vals, doy, year, kernels=kernels)

        # Calculate Saturn center to ring intercept vector
        if verbose:
            print('\tCalculating ring intercept event time and vector...')
        rho_vec_vals, t_ret_et_vals = cog.calc_rho_vec_km(t_oet_et_vals, planet,
                spacecraft, dsn, kernels=kernels)

        rho_km_vals = [spice.vnorm(vec) for vec in rho_vec_vals]



        # Interpolate to get sky frequency
        f_sky_hz_vals = np.interp(t_oet_spm_vals, f_spm, f_sky)


        # Calculate spacecraft event time
        if verbose:
            print('\tCalculating spacecraft event time...')
        t_set_et_vals = cog.calc_set_et(t_oet_et_vals, spacecraft, dsn)

        # Retrieve Saturn pole unit vector
        if verbose:
            print('\tRetrieving Saturn pole unit vector...')
        nhat_p = cog.get_pole(t_set_et_vals[0], planet)

        if verbose:
            print('\tCalculating ring longitude and ring azimuth...')

        # Calculated ring longitude and observed ring azimuth
        phi_rl_deg_vals, phi_ora_deg_vals = cog.calc_phi_deg(t_oet_et_vals,
                rho_vec_vals, spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating distance from spacecraft to ring intercept '
                + 'point...')
        # Calculate distance from spacecraft to ring intercept point
        D_km_vals = cog.calc_D_km(t_ret_et_vals, t_set_et_vals)

        if verbose:
            print('\tCalculating ring opening angle...')

        # Calculate ring opening angle
        B_deg_vals = cog.calc_B_deg(t_oet_et_vals, spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating Fresnel scale...')

        # Calculate Fresnel scale
        F_km_vals = cog.calc_F_km(D_km_vals, f_sky_hz_vals, B_deg_vals,
                phi_ora_deg_vals)

        t_ret_spm_vals = et_to_spm(t_ret_et_vals)
        t_set_spm_vals = et_to_spm(t_set_et_vals)

        if verbose:
            print('\tCalculating ring intercept velocities...')
        # Calculate ring intercept velocities
        rho_dot_kms_vals, phi_rl_dot_kms_vals = cog.calc_rip_velocity(
                rho_km_vals, phi_rl_deg_vals, step)

        if verbose:
            print('\tCalculating spacecraft state vector...')
        # Calculate spacecraft state vector
        R_sc_km_vals, R_sc_dot_kms_vals = cog.calc_sc_state(t_set_et_vals,
                spacecraft, planet, dsn, nhat_p)

        if verbose:
            print('\tCalculating impact radius...')
        # Calculate impact radius
        R_imp_km_vals = cog.calc_impact_radius_km(R_sc_km_vals, t_set_et_vals,
                spacecraft, dsn, nhat_p)

        if verbose:
            print('\tCalculating elevation angle...')

        # Calculate target angle above the horizon
        elev_deg_vals = cog.calc_elevation_deg(t_oet_et_vals, spacecraft, dsn)

        # Calculate beta
        if verbose:
            print('\tCalculating beta...')
        beta_vals = cog.calc_beta(B_deg_vals, phi_ora_deg_vals)
        B_eff_deg_vals = cog.calc_B_eff_deg(B_deg_vals, phi_ora_deg_vals)

        # Calculate when signal passes atmosphere + ionosphere
        ionos_occ_et_vals = cog.get_planet_occ_times(t_oet_et_vals, dsn,
                planet, spacecraft, height_above=5000.)
        self.ionos_occ_spm_vals = et_to_spm(ionos_occ_et_vals)
        self.ionos_occ_et_vals = ionos_occ_et_vals

        atmos_occ_et_vals = cog.get_planet_occ_times(t_oet_et_vals, dsn,
                planet, spacecraft, height_above=500.)
        self.atmos_occ_spm_vals = et_to_spm(atmos_occ_et_vals)


        # Set attributes, first block contains the inputs to *GEO.TAB
        self.t_oet_spm_vals = np.asarray(t_oet_spm_vals)
        self.t_ret_spm_vals = np.asarray(t_ret_spm_vals)
        self.t_set_spm_vals = np.asarray(t_set_spm_vals)

        self.t_set_et_vals = np.asarray(t_set_et_vals)
        self.t_oet_et_vals = np.asarray(t_oet_et_vals)

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
        self.elev_deg_vals = np.asarray(elev_deg_vals)
        self.naif_toolkit_version = self.__get_naif_version()
        self.beta_vals = np.asarray(beta_vals)
        self.B_eff_deg_vals = np.asarray(B_eff_deg_vals)

        # Add prof_dir entry to rev_info dict
        self.rev_info['prof_dir'] = self.get_profile_dir()

        if self.rev_info['prof_dir'] == '"BOTH"':
            self.split_ind = self.get_chord_ind()
        else:
            self.split_ind = None

        if len(self.atmos_occ_spm_vals) == 0:
            self.rev_info['planetary_occ_flag'] = '"N"'
        else:
            self.rev_info['planetary_occ_flag'] = '"Y"'



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


        self.freespace_km, self.freespace_spm = cog.get_freespace(
                t_ret_spm_vals, year, doy, rho_km_vals,
                phi_rl_deg_vals, t_oet_spm_vals, self.atmos_occ_spm_vals,
                split_ind=self.split_ind, kernels=kernels)

        # Write output data and label file if set
        if write_file:
            write_output_files(self)

    def __get_naif_version(self):
        """
        Return NAIF toolkit version used.

        Returns:
            :naif_str (*str*): NAIF toolkit version in format '"V.N0066"'.

        """
        naif_ver = spice.tkvrsn('TOOLKIT')
        naif_str = ('"V.' + naif_ver.split('_')[-1] + '"')
        return naif_str

    def get_profile_dir(self):
        """
        Purpose:
            Return observed profile direction.

        Returns:
            :prof_dir (*str*): Profile direction as '"INGRESS"',
                                 '"EGRESS"', or '"BOTH"'.
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
        """
        Purpose:
            Verify that an occultation with an increasing and decreasing
            radial velocity is actually a chord occultation.

        Returns:
            :prof_dir (*str*): Profile direction as
                                '"INGRESS"', '"EGRESS"', or '"BOTH"'.
        """
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

            # Remove false chord portions
            if self.verbose:
                print('\tRemoving portion blocked by atmosphere...')
            self.__remove_atmos_values()

        elif len(egr_blocked) == n_egr:
            prof_dir = '"INGRESS"'

            # Remove falsechord portions
            if self.verbose:
                print('\tRemoving portion blocked by atmosphere...')
            self.__remove_atmos_values()
        else:
            prof_dir = '"BOTH"'

        return prof_dir

    def __remove_atmos_values(self):
        npts = len(self.t_oet_spm_vals)

        # Index array of where atmosphere is blocking signal
        ind1 = np.argwhere(
                self.t_oet_spm_vals == min(self.atmos_occ_spm_vals))[0][0]
        ind2 = np.argwhere(
                self.t_oet_spm_vals == max(self.atmos_occ_spm_vals))[0][0]

        rind = np.arange(ind1, ind2+1, step=1)

        # Make sure ind2 > ind1
        if ind1 > ind2:
            ind1t = ind1
            ind2t = ind2t

            ind1 = ind2t
            ind2 = ind1t


        for attr, value in self.__dict__.items():
            if type(value) is not bool and len(value) == npts:
                setattr(self, attr, np.delete(value, rind))
            else:
                continue


    def get_chord_ind(self):
        """
        Purpose:
            Return index of where radial velocity sign change occurs in a chord
            occultation.

        Returns:
            :ind (*int*): Index of where chord occultation goes from '"INGRESS"'
                to '"EGRESS"' or vice versa.
        """

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

"""
Revisions:
"""
