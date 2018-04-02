#!/usr/bin/env python
"""

norm_diff_class.py

Purpose: Apply previously made calibration file to calibrate raw data
         (frequency correct and normalize), which produces the diffraction
         pattern given to Fresnel Inversion code.

Revisions:
      gjs_produce_normalized_diffraction_pattern.py
   Mar 07 2018 - gsteranka - Original version. Tentative rough draft, since
                             haven't yet chosen format of calibration file, and
                             also don't have Ryan's code to find window width
                             from a given window type and resolution.
   Mar 19 2018 - gsteranka - Edited to take new tentative final form of cal
                             file, which has the columns as 1) SPM 2) Predicted
                             sky frequency 3) Residual frequency fit
                             4) Freespace power spline fit 5) Frequency offset
                             fit. Also edited to save an obs file and return
                             everything in the obs file too
      produce_normalized_diffraction_pattern.py
   Mar 20 2018 - gsteranka - Eliminate debug steps, and switched obs_file to be
                             optional inputs
      norm_diff_class.py
   Mar 22 2018 - gsteranka - Edit previous code to be a class, with Fresnel
                            inversion inputs as attributes
   Mar 30 2018 - gsteranka - Edited reading in of geometry file to read in csv
                             format
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import spiceypy as spice
import sys

from rsr_reader import RSRReader
from resample_IQ import resample_IQ

class NormDiff(object):
    """Produce a normalized diffraction pattern using a previously made
    calibration file (cal_file). Contains everything needed for Fresnel
    Inversion as attributes.

    Attributes:
        rho_km_desired (np.ndarray): Ring intercept radius values, in km,
            at final spacing specified in dr_km
        spm_desired (np.ndarray): Set of SPM values corresponding to
            rho_km_desired
        p_norm_vals (np.ndarray): Power normalized to 1. This is the diffraction
            pattern that is input to the Fresnel inversion step
        phase_rad_vals (np.ndarray): Phase of the complex signal, in radians.
            This is the other part of the diffraction pattern that is input to
            the Fresnel Inversion step
        B_rad_vals (np.ndarray): Ring opening angle of Saturn
        D_km_vals (np.ndarray): Spacecraft-RIP (Ring Intercept Point) distance
        F_km_vals (np.ndarray): Fresnel scale, in km
        f_sky_pred_vals (np.ndarray): Predicted sky frequency, in Hz.
        phi_rad_vals (np.ndarray): Observed ring azimuth, in radians
        rho_dot_kms_vals (np.ndarray): Ring intercept point velocity
    """

    def __init__(self, rsr_file, dr_km, rho_km_range, res_km, window_type,
             geo_file, cal_file, dr_km_tol=0.01,
             decimate_16khz_to_1khz=False):
        """
        Instantiation defines all attributes of instance.

        Args:
            rsr_file (str): Full path name of RSR file to process
            dr_km (float): Desired final radial spacing of processed data, in km
            rho_km_range (list): Radius range to process over, in km
            res_km (float): Resolution to later do Fresnel Inversion at, in km
            window_type (str): Window type to later do Fresnel Inversion at
            geo_file (str): Full path name of a geoemtry file corresponding to
            rsr_file
            cal_file (str): Full path name of a calibration file corresponding to
                rsr_file
            dr_km_tol (float): Optional keyword argument, in km, that specifies the
                maximum distance the starting point of the final set of radius
                values can be away from an integer number of dr_km. For example, if
                you say dr_km_tol=0.01 and dr_km=0.25, the starting point of the
                final set of rho values could be 70000.26 km, but not 70000.261 km
            decimate_16khz_to_1khz (bool): Optional boolean keyword argument that
                specifies whether or not to decimate a 16khz sample rate file down
                to 1khz sample rate. Can only specify when rsr_file is 16khz. Value
                of this should match what you said in RSRfile.get_IQ method when you
                initially made the calibration file. Defualt is False
        """

        # Magical code to find radius range necessary for Fresnel Inversion.
        # For now just use the input radius range
        rho_km_range_inversion = rho_km_range

        # Read geometry file for radius-SPM conversion
        geo = pd.read_csv(geo_file, header=None)
        spm_geo = np.asarray(geo[0][:])
        rho_km_geo = np.asarray(geo[3][:])
        B_rad_geo = np.asarray(geo[6][:]) * spice.rpd()
        D_km_geo = np.asarray(geo[7][:])
        phi_rad_geo = np.asarray(geo[5][:]) * spice.rpd()
        rho_dot_kms_geo = np.asarray(geo[8][:])
        F_km_geo = np.asarray(geo[10][:])

        # Function to convert SPM to rho
        rho_interp_func = interp1d(spm_geo, rho_km_geo)

        # Find SPM range corresponding to rho range
        spm_start = max(spm_geo[rho_km_geo <= rho_km_range[0]])
        spm_end = min(spm_geo[rho_km_geo >= rho_km_range[1]])
        spm_range = [spm_start, spm_end]

        rsr_inst = RSRReader(rsr_file)
        (spm_vals, IQ_m) = rsr_inst.get_IQ(spm_range=spm_range,
                                           decimate_16khz_to_1khz=
                                           decimate_16khz_to_1khz)
        rho_km_vals = rho_interp_func(spm_vals)

        cal = pd.read_csv(cal_file, delim_whitespace=True, header=None)

        spm_cal = np.asarray(cal[0][:])
        f_sky_pred_cal = np.asarray(cal[1][:])
        p_free_cal = np.asarray(cal[3][:])
        f_offset_fit_cal = np.asarray(cal[4][:])
        rho_km_cal = rho_interp_func(spm_cal)

        # If specified rho range that cal file can't cover
        if (spm_range[0]<min(spm_cal)) | (spm_range[1]>max(spm_cal)):
            print('ERROR (NormDiff): Specified cal file is missing \n'
                  +'points required to pre-process points in rho_km_range.\n'
                  +'Calibration file rho range is '+str(min(rho_km_cal))+'km, '
                  +str(max(rho_km_cal))+'km, while specified rho_km_range is \n'
                  +str(rho_km_range[0])+'km, '+str(rho_km_range[1])+'km')
            sys.exit()

        # Inteprolate frequency offset to finer spacing in preparation
        # for integration
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

        (rho_km_desired, IQ_c_desired) = resample_IQ(rho_km_vals, IQ_c, dr_km,
                                                     dr_km_tol=dr_km_tol)

        (rho_km_desired, IQ_m_desired) = resample_IQ(rho_km_vals, IQ_m, dr_km,
                                                     dr_km_tol=dr_km_tol)

        # Interpolate freespace power (spline fit) to rho_km_desired
        p_free_interp_func = interp1d(rho_km_cal, p_free_cal)
        p_free = p_free_interp_func(rho_km_desired)

        p_norm_vals = (abs(IQ_m_desired)**2)/p_free
        phase_rad_vals = np.arctan2(np.imag(IQ_c_desired),
                                    np.real(IQ_c_desired))

        spm_desired_func = interp1d(rho_km_geo, spm_geo)
        spm_desired = spm_desired_func(rho_km_desired)

        self.__set_attributes(rho_km_desired, spm_desired, p_norm_vals,
                              phase_rad_vals,
                              rho_km_geo, B_rad_geo, D_km_geo, F_km_geo,
                              phi_rad_geo, rho_dot_kms_geo,
                              rho_km_cal, f_sky_pred_cal)

    def __set_attributes(self, rho_km_desired, spm_desired, p_norm_vals,
                         phase_rad_vals,
                         rho_km_geo, B_rad_geo, D_km_geo, F_km_geo, phi_rad_geo,
                         rho_dot_kms_geo,
                         rho_km_cal, f_sky_pred_cal):
        """Private method called by __init__ to set attributes of the
        instance"""

        self.rho_km_desired = rho_km_desired
        self.spm_desired = spm_desired
        self.p_norm_vals = p_norm_vals
        self.phase_rad_vals = phase_rad_vals

        # Ring opening angle at final spacing
        B_rad_func = interp1d(rho_km_geo, B_rad_geo, fill_value='extrapolate')
        self.B_rad_vals = B_rad_func(rho_km_desired)

        # Spacecraft - RIP distance at final spacing
        D_km_func = interp1d(rho_km_geo, D_km_geo, fill_value='extrapolate')
        self.D_km_vals = D_km_func(rho_km_desired)

        # Fresnel scale at final spacing
        F_km_func = interp1d(rho_km_geo, F_km_geo, fill_value='extrapolate')
        self.F_km_vals = F_km_func(rho_km_desired)

        # Sky frequency at final spacing
        f_sky_pred_func = interp1d(rho_km_cal, f_sky_pred_cal)
        self.f_sky_pred_vals = f_sky_pred_func(rho_km_desired)

        # Observed Ring Azimuth at final spacing
        phi_rad_func = interp1d(rho_km_geo, phi_rad_geo, fill_value='extrapolate')
        self.phi_rad_vals = phi_rad_func(rho_km_desired)

        # RIP velocity at final spacing
        rho_dot_kms_func = interp1d(rho_km_geo, rho_dot_kms_geo,
                                    fill_value='extrapolate')
        self.rho_dot_kms_vals = rho_dot_kms_func(rho_km_desired)


    def save_obs(self, obs_file):
        """Save an obs file, which contains all of the attributes in the
        instance, and therefore has everything needed for Fresnel Inversion

        Args:
            obs_file (str): Full path name of obs file to be created
        """

        np.savetxt(obs_file,
                   np.c_[self.rho_km_desired, self.spm_desired,
                         self.p_norm_vals, self.phase_rad_vals,
                         self.B_rad_vals, self.D_km_vals, self.F_km_vals,
                         self.f_sky_pred_vals, self.phi_rad_vals,
                         self.rho_dot_kms_vals],
                   fmt='%32.16f '*10)


if __name__ == '__main__':
    pass
