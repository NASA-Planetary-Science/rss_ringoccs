#!/usr/bin/env python
"""

norm_diff_class.py

Purpose: Apply previously made calibration file to calibrate raw data
         (frequency correct and normalize), which produces the diffraction
         pattern given to Fresnel Inversion code.

Revisions:
      gjs_produce_normalized_diffraction_pattern.py
   2018 Mar 07 - gsteranka - Original version. Tentative rough draft, since
                             haven't yet chosen format of calibration file, and
                             also don't have Ryan's code to find window width
                             from a given window type and resolution.
   2018 Mar 19 - gsteranka - Edited to take new tentative final form of cal
                             file, which has the columns as 1) SPM 2) Predicted
                             sky frequency 3) Residual frequency fit
                             4) Freespace power spline fit 5) Frequency offset
                             fit. Also edited to save an obs file and return
                             everything in the obs file too
      produce_normalized_diffraction_pattern.py
   2018 Mar 20 - gsteranka - Eliminate debug steps, and switched obs_file to be
                             optional inputs
      norm_diff_class.py
   2018 Mar 22 - gsteranka - Edit previous code to be a class, with Fresnel
                            inversion inputs as attributes
   2018 Mar 30 - gsteranka - Edited reading in of geometry file to read in csv
                             format
   2018 Apr 03 - gsteranka - When first defining p_norm_vals, do it by
                             normalizing IQ_c_desired, not IQ_m_desired
   2018 Apr 09 - gsteranka - Obs files now save absolute value of rho dot,
                             which eliminates a lot of error checking in
                             Fresnel Inversion step. Added
                             "fill_value='extrapolate'" to rho_interp_func,
                             p_free_interp_func, spm_desired_func, and
                             f_sky_pred_func
   2018 Apr 19 - gsteranka - Changed attribute names to agree w/ what Ryan
                             uses: rho_km_desired --> rho_km_vals
                                   spm_desired --> spm_vals
                                   f_sky_pred_vals --> f_sky_hz_vals
   2018 Apr 30 - gsteranka - Edited several portions to be compatible with
                             chord occultations. Changes include:
                                 1) Add "is_chord" keyword to __init__
                                 2) Set 1st element of rho_dot_kms_geo equal
                                    to 2nd element, since 1st element is off
                                    due to how it's originally calculated
                                 3) Add if statement separating chord/non-chord
                                    cases
                                 4) __set_attributes() metho takes spm_geo
                                 5) Added "end_of_chord_ing" attribute
                                 6) save_obs() method uses "end_of_chord_ing"
                                    attribute to separate ingress/egress
                                    portions of chord occultation and make
                                    separate files
                                 7) Separate case for chord occultations to
                                    convert "rho_km_range" to "spm_range"
   2018 May 07 - gsteranka - Edited how it detects where ingress switches to
                             egress in a chord occultation. It used to use
                             rho_dot, but I find this is liable to slight
                             interpolation error which causes it to be off by
                             a bit. I now use rho_km_vals directly
   2018 May 23 - gsteranka - Edited to read cal file as a csv, rather than
                             white-space delimited. Also added the
                             t_ret_et_vals, t_set_et_vals, and phi_rl_rad_vals
                             attributes
   2018 Jun 04 - gsteranka - Added history attribute and __set_history() method
   2018 Jun 11 - gsteranka - Accept rsr_inst, geo_inst, and cal_inst instead of
                             the files for all of those. "cal_inst" can be made
                             from Calibration class calibration_class.py if you
                             haven't made a cal file yet, or from MakeCalInst
                             in make_cal_inst.py if you have made a cal file
                             already
   2018 Jun 18 - gsteranka - Added kernels input. Only needs to be leap second
                             kernel
   2018 Jun 21 - gsteranka - Change spacecraft event time and ring event time
                             attributes to SPM instead of et. Removed kernels
                             input
   2018 Jun 27 - gsteranka - Adjust CAL file use so sky frequency column has
                             frequency offset fit added to it
   2018 Jun 28 - gsteranka - Change spm_vals attribute to t_oet_spm_vals
"""

import copy
import numpy as np
import os
import pdb
import platform
from scipy.interpolate import interp1d
import spiceypy as spice
import sys
import time

try:
    from rsr_reader import RSRReader
    from resample_IQ import resample_IQ
except ImportError:
    from ..rsr_reader.rsr_reader import RSRReader
    from .resample_IQ import resample_IQ


class NormDiff(object):
    """
    Purpose:
    Produce a normalized diffraction pattern using a previously made
    calibration file (cal_file). Contains everything needed for Fresnel
    Inversion as attributes.

    Attributes:
        rho_km_vals (np.ndarray): Ring intercept radius values, in km,
            at final spacing specified in dr_km
        t_oet_spm_vals (np.ndarray): Set of SPM values corresponding to
            rho_km_desired
        p_norm_vals (np.ndarray): Power normalized to 1. This is the diffraction
            pattern that is input to the Fresnel inversion step
        phase_rad_vals (np.ndarray): Phase of the complex signal, in radians.
            This is the other part of the diffraction pattern that is input to
            the Fresnel Inversion step
        B_rad_vals (np.ndarray): Ring opening angle of Saturn
        D_km_vals (np.ndarray): Spacecraft-RIP (Ring Intercept Point) distance
        F_km_vals (np.ndarray): Fresnel scale, in km
        f_sky_hz_vals (np.ndarray): Predicted sky frequency, in Hz.
        phi_rad_vals (np.ndarray): Observed ring azimuth, in radians
        phi_rl_rad_vals (np.ndarray): Ring longitude, in radians
        t_ret_spm_vals (np.ndarray): Ring event time in SPM
        t_set_spm_vals (np.ndarray): Spacecraft event time in SPM
        rho_dot_kms_vals (np.ndarray): Ring intercept point velocity
        rho_corr_pole_km_vals (np.ndarray): Radius correction due to Saturn's
            pole direction
        rho_corr_timing_km_vals (np.ndarray): Radius correction due to timing
            offset
        tau_threshold_vals (np.ndarray): Threshold optical depth assuming 1km
            desired resolution
        end_of_chord_ing (int): Index number of final ingress portion of chord
            occultation. Set to "None" if "is_chord" keyword is False
        history (dict): Recorded information about the run
    """

    def __init__(self, rsr_inst, dr_km, geo_inst, cal_inst,
            dr_km_tol=0.01, is_chord=False):
        """
        Purpose:
        Instantiation defines all attributes of instance.

        Args:
            rsr_inst: Instance of the RSRReader class. Linked to the full
                path name of RSR file to process
            dr_km (float): Desired final radial spacing of processed data, in km
            geo_inst (str): Instance of Geometry class linked to rsr_inst input
            cal_inst: Calibration instance linked to rsr_inst input. Made using
                Calibration class if you haven't made a cal_file yet, and made
                using MakeCalInst class if you have made a cal_file
            dr_km_tol (float): Optional keyword argument, in km, that specifies the
                maximum distance the starting point of the final set of radius
                values can be away from an integer number of dr_km. For example, if
                you say dr_km_tol=0.01 and dr_km=0.25, the starting point of the
                final set of rho values could be 70000.26 km, but not 70000.261 km
            is_chord (bool): Set as True if the occultation is a chord occultation

        Dependencies:
            [1] RSRReader
            [2] resample_IQ
            [3] spm_to_et
            [4] numpy
            [5] os
            [6] platform
            [7] scipy.interpolate
            [8] spiceypy
            [9] sys
            [10] time

        Warnings:
            [1] Be sure to set decimate_16khz_to_1khz=True when making rsr_inst
                if you did so before doing the previous processing steps to make
                GEO and CAL files. Otherwise, (1) it will take forever, and
                (2) the power will not normalize correctly.
            [2] If it's a chord occultation, you will get problems if you don't
                set is_chord=True
            [3] Radius correction due to pole direction and timing offset not
                yet implemented, so the corresponding attributes
                (rho_corr_pole_km_vals and rho_corr_timing_km_vals) are set
                to arrays of -999 right now
        """

        spm_geo = geo_inst.t_oet_spm_vals
        rho_km_geo = geo_inst.rho_km_vals
        rho_dot_kms_geo = geo_inst.rho_dot_kms_vals

        # First element has crazy value due to how rho_dot is calculated, so
        #     eliminate it in favor of next element
        rho_dot_kms_geo[0] = rho_dot_kms_geo[1]

        # Check if rho_dot has both positive and negative elements. If so, it
        #     might be a chord occultation, so let the user know if they didn't
        #     specify it as one. Ignore very first and last elements because
        #     those sometimes go to crazy values
        rho_dot_trim = rho_dot_kms_geo[1:-2]
        if (len(rho_dot_trim[rho_dot_trim > 0]) != 0) & (
                (len(rho_dot_trim[rho_dot_trim < 0]) != 0)) & (
                is_chord is False):
            print('WARNING: possible chord occultation detected, but is_chord '
                + 'is False. Did you forget to set is_chord=True?')

        # Function to convert SPM to rho
        rho_interp_func = interp1d(spm_geo, rho_km_geo, fill_value='extrapolate')

        # SPM range of geo_inst (pretty much always the full occultation)
        spm_range = [min(spm_geo), max(spm_geo)]

        spm_vals = rsr_inst.spm_vals
        IQ_m = rsr_inst.IQ_m
        rho_km_vals = rho_interp_func(spm_vals)

        spm_cal = cal_inst.t_oet_spm_vals
        f_sky_pred_cal = cal_inst.f_sky_hz_vals
        p_free_cal = cal_inst.p_free_vals
        dummy_spm, f_sky_pred_file = rsr_inst.get_f_sky_pred(f_spm=spm_cal)
        f_offset_fit_cal = f_sky_pred_cal - f_sky_pred_file
        rho_km_cal = rho_interp_func(spm_cal)

        # If that cal file can't cover entire SPM range from geo_inst
        if (spm_range[0]<min(spm_cal)) | (spm_range[1]>max(spm_cal)):
            print('\nWARNING (NormDiff): Specified cal file is missing '
                + 'points required to pre-process points in SPM range of '
                + 'geo_inst. Calibration SPM range is ['
                + str(min(spm_cal))+', ' + str(max(spm_cal))
                +'], while geo_inst SPM range is ['
                + str(min(spm_geo))+', '+str(max(spm_geo))+']. Reducing to '
                + 'SPM range of CAL file')
            _ind = ((spm_vals >= min(spm_cal)) &
                (spm_vals <= max(spm_cal))).nonzero()
            spm_vals = spm_vals[_ind]
            IQ_m = IQ_m[_ind]
            rho_km_vals = rho_km_vals[_ind]

        # Inteprolate frequency offset to finer spacing in preparation
        #     for integration
        dt = 0.1
        n_pts = round((spm_cal[-1] - spm_cal[0])/dt)
        spm_cal_interp = spm_cal[0] + dt*np.arange(n_pts)
        f_offset_fit_func = interp1d(spm_cal, f_offset_fit_cal,
            fill_value='extrapolate')
        f_offset_fit_interp = f_offset_fit_func(spm_cal_interp)

        # Integrate frequency offset to detrend IQ_m
        f_detrend_interp = np.cumsum(f_offset_fit_interp)*dt
        f_detrend_interp_rad = f_detrend_interp * spice.twopi()
        f_detrend_rad_func = interp1d(spm_cal_interp, f_detrend_interp_rad,
            fill_value='extrapolate')
        f_detrend_rad = f_detrend_rad_func(spm_vals)

        # Detrend raw measured I and Q
        IQ_c = IQ_m*np.exp(-1j*f_detrend_rad)

        if not is_chord:
            rho_km_desired, IQ_c_desired = resample_IQ(rho_km_vals, IQ_c, dr_km,
                dr_km_tol=dr_km_tol)

            # Interpolate freespace power (spline fit) to rho_km_desired
            p_free_interp_func = interp1d(rho_km_cal, p_free_cal,
                fill_value='extrapolate')
            p_free = p_free_interp_func(rho_km_desired)
            
            # Construct power from IQ_c_desired and not IQ_m_desired or
            #     resampled power, because p_free was constructed from
            #     downsampled IQ_c. It's important to not use IQ_m_desired
            #     because I and Q each vary sinusoidally, meaning that if you
            #     resample too much, they get averaged to near zero.
            p_norm_vals = (abs(IQ_c_desired)**2)/p_free
            phase_rad_vals = np.arctan2(np.imag(IQ_c_desired),
                np.real(IQ_c_desired))
            
            spm_desired_func = interp1d(rho_km_geo, spm_geo,
                fill_value='extrapolate')
            spm_desired = spm_desired_func(rho_km_desired)

            end_of_chord_ing = None

        else:

            rho_km_vals_diff = np.zeros(len(rho_km_vals))
            rho_km_vals_diff[1:] = np.diff(rho_km_vals)
            rho_km_vals_diff[0] = rho_km_vals_diff[1]
            ing_ind = rho_km_vals_diff < 0
            egr_ind = rho_km_vals_diff > 0

            rho_km_ing, IQ_c_ing = resample_IQ(rho_km_vals[ing_ind],
                IQ_c[ing_ind], dr_km, dr_km_tol=dr_km_tol)
            rho_km_egr, IQ_c_egr = resample_IQ(rho_km_vals[egr_ind],
                IQ_c[egr_ind], dr_km, dr_km_tol=dr_km_tol)

            spm_ing_func = interp1d(rho_km_vals[ing_ind], spm_vals[ing_ind],
                fill_value='extrapolate')
            spm_egr_func = interp1d(rho_km_vals[egr_ind], spm_vals[egr_ind],
                fill_value='extrapolate')
            spm_ing = spm_ing_func(rho_km_ing)
            spm_egr = spm_egr_func(rho_km_egr)

            p_free_ing = np.interp(spm_ing, spm_cal, p_free_cal)
            p_free_egr = np.interp(spm_egr, spm_cal, p_free_cal)

            p_norm_ing = (abs(IQ_c_ing)**2)/p_free_ing
            p_norm_egr = (abs(IQ_c_egr)**2)/p_free_egr
            phase_rad_ing = np.arctan2(np.imag(IQ_c_ing), np.real(IQ_c_ing))
            phase_rad_egr = np.arctan2(np.imag(IQ_c_egr), np.real(IQ_c_egr))

            rho_km_desired = np.concatenate((rho_km_ing, rho_km_egr))
            spm_desired = np.concatenate((spm_ing, spm_egr))
            p_norm_vals = np.concatenate((p_norm_ing, p_norm_egr))
            phase_rad_vals = np.concatenate((phase_rad_ing, phase_rad_egr))

            end_of_chord_ing = len(rho_km_ing) - 1


        self.__set_history(rsr_inst, dr_km, geo_inst, cal_inst, dr_km_tol)
        self.__set_attributes(rho_km_desired, spm_desired, p_norm_vals,
            phase_rad_vals,
            spm_geo, rho_dot_kms_geo, geo_inst,
            spm_cal, f_sky_pred_cal, rsr_inst,
            end_of_chord_ing=end_of_chord_ing)


    def __set_attributes(self, rho_km_desired, spm_desired, p_norm_vals,
            phase_rad_vals,
            spm_geo, rho_dot_kms_geo, geo_inst,
            spm_cal, f_sky_pred_cal, rsr_inst,
            end_of_chord_ing=None):
        """
        Private method called by __init__ to set attributes of the
        instance
        """

        self.rho_km_vals = rho_km_desired
        self.t_oet_spm_vals = spm_desired
        self.p_norm_vals = p_norm_vals
        self.phase_rad_vals = phase_rad_vals

        B_rad_geo = geo_inst.B_deg_vals*spice.rpd()
        D_km_geo = geo_inst.D_km_vals
        phi_rad_geo = geo_inst.phi_ora_deg_vals*spice.rpd()
        F_km_geo = geo_inst.F_km_vals
        t_ret_geo = geo_inst.t_ret_spm_vals
        t_set_geo = geo_inst.t_set_spm_vals
        phi_rl_rad_geo = geo_inst.phi_rl_deg_vals*spice.rpd()
        #rho_corr_pole_km_geo = geo_inst.rho_corr_pole_km_vals
        #rho_corr_timing_km_geo = geo_inst.rho_corr_timing_km_vals

        # Ring opening angle at final spacing
        B_rad_func = interp1d(spm_geo, B_rad_geo, fill_value='extrapolate')
        self.B_rad_vals = B_rad_func(spm_desired)

        # Spacecraft - RIP distance at final spacing
        D_km_func = interp1d(spm_geo, D_km_geo, fill_value='extrapolate')
        self.D_km_vals = D_km_func(spm_desired)

        # Fresnel scale at final spacing
        F_km_func = interp1d(spm_geo, F_km_geo, fill_value='extrapolate')
        self.F_km_vals = F_km_func(spm_desired)

        # Sky frequency at final spacing
        f_sky_pred_func = interp1d(spm_cal, f_sky_pred_cal,
            fill_value='extrapolate')
        self.f_sky_hz_vals = f_sky_pred_func(spm_desired)

        # Observed Ring Azimuth at final spacing
        phi_rad_func = interp1d(spm_geo, phi_rad_geo, fill_value='extrapolate')
        self.phi_rad_vals = phi_rad_func(spm_desired)

        # Ring Event Time at final spacing
        t_ret_func = interp1d(spm_geo, t_ret_geo, fill_value='extrapolate')
        self.t_ret_spm_vals = t_ret_func(spm_desired)

        # Spacecraft Event Time at final spacing
        t_set_func = interp1d(spm_geo, t_set_geo, fill_value='extrapolate')
        self.t_set_spm_vals = t_set_func(spm_desired)

        # Ring longitude at final spacing
        phi_rl_rad_func = interp1d(spm_geo, phi_rl_rad_geo,
            fill_value='extrapolate')
        self.phi_rl_rad_vals = phi_rl_rad_func(spm_desired)

        # RIP velocity at final spacing
        rho_dot_kms_func = interp1d(spm_geo, rho_dot_kms_geo,
            fill_value='extrapolate')
        self.rho_dot_kms_vals = rho_dot_kms_func(spm_desired)

        # Radius correction due to Saturn's pole direction at final spacing
        #rho_corr_pole_km_func = interp1d(spm_geo, rho_corr_pole_km_geo,
        #    fill_value='extrapolate')
        #self.rho_corr_pole_km_vals = rho_corr_pole_km_func(spm_desired)
        self.rho_corr_pole_km_vals = np.zeros(len(spm_desired))

        # Radius correction due to timing offset at final spacing
        #rho_corr_timing_km_func = interp1d(spm_geo, rho_corr_timing_km_geo,
        #    fill_value='extrapolate')
        #self.rho_corr_timing_km_vals = rho_corr_timing_km_func(spm_desired)
        self.rho_corr_timing_km_vals = np.zeros(len(spm_desired))

        self.tau_threshold_vals = np.zeros(len(spm_desired))

        self.end_of_chord_ing = end_of_chord_ing


    def __set_history(self, rsr_inst, dr_km, geo_inst, cal_inst, dr_km_tol):

        input_var_dict = {'rsr_inst': rsr_inst.history, 'dr_km': dr_km,
            'geo_file': geo_inst.history, 'cal_inst': cal_inst.history}
        input_kw_dict = {'dr_km_tol': dr_km_tol}
        hist_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__.split('/')[-1],
            'Source Directory': __file__.rsplit('/',1)[0] +'/',
            'Input Variables': input_var_dict,
            'Input Keywords':input_kw_dict}
        self.history = hist_dict


    def chord_split(self):
        """
        Split chord occultation instance into an ingress instance and an egress
        instance, and return them.

        Returns:
            norm_diff_inst_ing: Ingress instance
            norm_diff_inst_egr: Egress instance

        Notes:
            [1] Is it bad practice to create separate copies of the instances
                and return them? Not sure what else to do.
        """

        if self.end_of_chord_ing is None:
            print('WARNING (NormDiff.chord_split()): Not a chord occultation, '+
                'so can\'t split. Ignoring call to this method, returning None')
            return None

        # Record all attribute names
        rho_km_vals = self.rho_km_vals
        t_oet_spm_vals = self.t_oet_spm_vals
        p_norm_vals = self.p_norm_vals
        phase_rad_vals = self.phase_rad_vals
        B_rad_vals = self.B_rad_vals
        D_km_vals = self.D_km_vals
        F_km_vals = self.F_km_vals
        f_sky_hz_vals = self.f_sky_hz_vals
        phi_rad_vals = self.phi_rad_vals
        phi_rl_rad_vals = self.phi_rl_rad_vals
        t_ret_spm_vals = self.t_ret_spm_vals
        t_set_spm_vals = self.t_set_spm_vals
        rho_dot_kms_vals = self.rho_dot_kms_vals
        rho_corr_pole_km_vals = self.rho_corr_pole_km_vals
        rho_corr_timing_km_vals = self.rho_corr_timing_km_vals
        tau_threshold_vals = self.tau_threshold_vals
        end_of_chord_ing = self.end_of_chord_ing
        history = self.history

        norm_diff_inst_ing = copy.deepcopy(self)
        norm_diff_inst_egr = copy.deepcopy(self)

        ind = end_of_chord_ing

        # Ingress instance
        norm_diff_inst_ing.rho_km_vals = rho_km_vals[0:ind]
        norm_diff_inst_ing.t_oet_spm_vals = (
            t_oet_spm_vals[0:ind])
        norm_diff_inst_ing.p_norm_vals = p_norm_vals[0:ind]
        norm_diff_inst_ing.phase_rad_vals = phase_rad_vals[0:ind]
        norm_diff_inst_ing.B_rad_vals = B_rad_vals[0:ind]
        norm_diff_inst_ing.D_km_vals = D_km_vals[0:ind]
        norm_diff_inst_ing.F_km_vals = F_km_vals[0:ind]
        norm_diff_inst_ing.f_sky_hz_vals = (
            f_sky_hz_vals[0:ind])
        norm_diff_inst_ing.phi_rad_vals = phi_rad_vals[0:ind]
        norm_diff_inst_ing.phi_rl_rad_vals = (
            phi_rl_rad_vals[0:ind])
        norm_diff_inst_ing.t_ret_spm_vals = (
            t_ret_spm_vals[0:ind])
        norm_diff_inst_ing.t_set_spm_vals = (
            t_set_spm_vals[0:ind])
        norm_diff_inst_ing.rho_dot_kms_vals = (
            rho_dot_kms_vals[0:ind])
        norm_diff_inst_ing.rho_corr_pole_km_vals = (
            rho_corr_pole_km_vals[0:ind])
        norm_diff_inst_ing.rho_corr_timing_km_vals = (
            rho_corr_timing_km_vals[0:ind])
        norm_diff_inst_ing.tau_threshold_vals = (
            tau_threshold_vals[0:ind])
        norm_diff_inst_ing.end_of_chord_ing = None
        norm_diff_inst_ing.history = history

        # Egress instance
        norm_diff_inst_egr.rho_km_vals = rho_km_vals[ind+1:]
        norm_diff_inst_egr.t_oet_spm_vals = (
            t_oet_spm_vals[ind+1:])
        norm_diff_inst_egr.p_norm_vals = p_norm_vals[ind+1:]
        norm_diff_inst_egr.phase_rad_vals = phase_rad_vals[ind+1:]
        norm_diff_inst_egr.B_rad_vals = B_rad_vals[ind+1:]
        norm_diff_inst_egr.D_km_vals = D_km_vals[ind+1:]
        norm_diff_inst_egr.F_km_vals = F_km_vals[ind+1:]
        norm_diff_inst_egr.f_sky_hz_vals = (
            f_sky_hz_vals[ind+1:])
        norm_diff_inst_egr.phi_rad_vals = phi_rad_vals[ind+1:]
        norm_diff_inst_egr.phi_rl_rad_vals = (
            phi_rl_rad_vals[ind+1:])
        norm_diff_inst_egr.t_ret_spm_vals = (
            t_ret_spm_vals[ind+1:])
        norm_diff_inst_egr.t_set_spm_vals = (
            t_set_spm_vals[ind+1:])
        norm_diff_inst_egr.rho_dot_kms_vals = (
            rho_dot_kms_vals[ind+1:])
        norm_diff_inst_egr.rho_corr_pole_km_vals = (
            rho_corr_pole_km_vals[ind+1:])
        norm_diff_inst_egr.rho_corr_timing_km_vals = (
            rho_corr_timing_km_vals[ind+1:])
        norm_diff_inst_egr.tau_threshold_vals = (
            tau_threshold_vals[ind+1:])
        norm_diff_inst_egr.end_of_chord_ing = None
        norm_diff_inst_egr.history = history

        return norm_diff_inst_ing, norm_diff_inst_egr


if __name__ == '__main__':
    pass
