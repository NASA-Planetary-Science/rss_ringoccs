#!/usr/bin/env python
"""

freq_offset_fit.py

Purpose: Makes a fit to the frequency offset made from freq_offset.py, using
         the frequency offset, predicted sky frequency, reconstructed sky
         frequency, and a fit to residual frequency.

Revisions:
      gjs_freq_offset_fit.py
   2018 Feb 21 - gsteranka - Original version
      gjs_freq_offset_fit_v2.py
   2018 Mar 06 - gsteranka - Edited to define all attributes during
                             instantiation
   2018 Mar 07 - gsteranka - Eliminated freq_offset_inst input in favor of
                             f_SPM and f_offset inputs. This is consistent
                             with how gjs_freq_offset_v2.py now works
      freq_offset_fit.py
   2018 Mar 20 - gsteranka - Copy to official version and remove debug steps
   2018 Mar 20 - gsteranka - Added get_f_offset_fit method, and an error check
                             if all points fall inside of rho exclusion zone
   2018 Apr 27 - gsteranka - Added TEST keyword to __init__ and
                             set_f_sky_resid_fit methods to plot observed
                             residual frequency, points used as freespace to
                             fit over, and the resulting fit
   2018 May 03 - gsteranka - Check if behind Saturn's ionosphere before fitting
   2018 May 08 - gsteranka - Added "spm_include" keyword to set_f_sky_resid_fit
                             method to optionally replace the rho_exclude
                             input. Added GUI for residual frequency fit, and
                             it defaults to using it
   2018 May 30 - gsteranka - Added _spm_include and _k attributes, and edited
                             use of GUI to make initial guess at fit. Eliminated
                             rho_exclude keyword
   2018 Jun 05 - gsteranka - Added sc_name keyword. Default is Cassini
"""

import numpy as np
from numpy.polynomial import polynomial as poly
import os
import platform
from scipy.interpolate import interp1d
import spiceypy.spiceypy as spice
import sys
import time

from tkinter import Tk

try:
    from calc_f_sky_recon import calc_f_sky_recon
    from cassini_blocked import cassini_blocked
    from f_resid_fit_gui import FResidFitGui
except ImportError:
    from .calc_f_sky_recon import calc_f_sky_recon
    from ..tools.cassini_blocked import cassini_blocked
    from .f_resid_fit_gui import FResidFitGui

class FreqOffsetFit(object):
    """Class to make a fit to extracted frequency offset. Uses predicted sky
    frequency, reconstructed sky frequency, and a fit to residual sky frequency
    to do so.

    Example:
        >>> rsr_inst = RSRReader(rsr_file)
        >>> (f_spm, f_offset) = (f_spm, f_offset) = \
                calc_freq_offset(spm_vals, IQ_m)
        >>> fit_inst = FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset, \
                f_uso, kernels, k=k, rho_exclude=rho_exclude)
        >>> (spm_vals, IQ_c) = fit_inst.get_IQ_c(spm_vals=spm_vals, IQ_m=IQ_m)

    Attributes:
        _spm_include (list): Linked to the input range spm_include
        _k (int): Linked to the input fit order k
        __geo_inst: Instance of Geometry class
        __f_uso (float): USO frequency used
        __USE_GUI (bool): Linked to USE_GUI input
        __f_spm (np.ndarray): SPM values that all sky frequencies are
            evaluated at
        __f_rho (np.ndarray): Rho values that all sky frequencies are
            evaluated at
        __f_offset (np.ndarray): Extracted frequency offset
        __f_offset_fit (np.ndarray): Fit to extracted frequency offset
        __f_sky_pred (np.ndarray): Predicted sky frequency
        __f_sky_recon (np.ndarray): Reconstructed sky frequency
        __f_sky_resid_fit (np.ndarray): Fit to observed residual frequency
        __spm_vals (np.ndarray): SPM values at raw resolution over
            full data set
        __IQ_m (np.ndarray): Raw measured complex signal over full data set
        history (dict): Recorded information about the run
        """

    def __init__(self, rsr_inst, geo_inst, f_spm, f_offset,
            f_uso, k=9, spm_include=None, sc_name='Cassini',
            USE_GUI=True, TEST=False):
        """
        Purpose:
        Define all attributes associated with data set, and make a fit to
        frequency offset using the default parameters

        Args:
            rsr_inst: Instance of the RSRReader class
            geo_inst: Instance of geometry class. Contains attributes t_oet_spm
                and rho_km_vals
            f_spm (np.ndarray): SPM array that frequency offset was
                extracted at
            f_offset (np.ndarray): Extracted frequency offset
            f_uso (float): USO frequency for the data set
            k (int): Order of the polynomial fit made to residual frequency
            rho_exclude (list): Set of radius regions to exclude when making
                fit to residual frequency. Specify in km. Default is to
                exclude B ring region
            sc_name (str): Name of spacecraft. Default is 'Cassini'
            spm_include (list): Set of SPM regions to include when making fit
                to residual frequency. By default, only rho_exclude is used.
                If this keyword is specified, it overrides anything input to
                rho_exclude keyword. This is meant as an optional replacement
                of rho_exclude
            USE_GUI (bool): Use the interactive GUI to make a residual
                frequency fit. This is highly recommended
            TEST (bool): Optional test plot

        Dependencies:
            [1] RSRReader
            [2] Geometry
            [3] numpy
            [4] os
            [5] platform
            [6] scipy.interpolate
            [7] sys
            [8] time

        Warnings:
            [1] If you don't use the GUI the first time you run a data set,
                then the resulting residual frequency fit is liable to be bad
                due to eccentric ringlets
            [2] Inputting an incorrect or incomplete set of kernels to Geometry
                routine before running this will give you a cryptic error
                message from spiceypy routines
        """

        # Keep the RSRReader instance, kernels, Geometry instance, and f_USO
        #     as attributes
        self.__rsr_inst = rsr_inst
        self.__kernels = geo_inst.kernels
        self.__geo_inst = geo_inst
        self.__f_uso = f_uso

        # Attributes for components of frequency offset, except for residual
        # frequency and its fit
        self.__f_spm = np.asarray(f_spm)
        self.__f_offset = np.asarray(f_offset)
        f_spm, self.__f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=self.__f_spm)
        self.__f_sky_recon = calc_f_sky_recon(self.__f_spm, rsr_inst, sc_name,
            f_uso, self.__kernels)

        # Interpolate geometry file rho's to rho's for f_spm and
        # spm_vals (raw resolution)
        spm_geo = np.asarray(geo_inst.t_oet_spm_vals)
        rho_geo = np.asarray(geo_inst.rho_km_vals)
        rho_interp_func = interp1d(spm_geo, rho_geo, fill_value='extrapolate')
        self.__f_rho = rho_interp_func(self.__f_spm)

        # Attribute keeping track of fit parameter
        self._spm_include = spm_include

        # Set attributes for residual frequency fit, and the
        # new frequency offset fit
#         self.set_f_sky_resid_fit(k=k, rho_exclude=rho_exclude,
#             spm_include=spm_include, USE_GUI=USE_GUI, TEST=TEST)
        self.set_f_sky_resid_fit(k=k, spm_include=spm_include, USE_GUI=USE_GUI,
            TEST=TEST)

        # Get I and Q from RSR object
        self.__spm_vals = rsr_inst.spm_vals
        self.__IQ_m = rsr_inst.IQ_m


#     def set_f_sky_resid_fit(self, k=9, rho_exclude=None, spm_include=None,
#             USE_GUI=True, TEST=False)
    def set_f_sky_resid_fit(self, k=9, spm_include=None, USE_GUI=True,
            TEST=False):
        """
        Purpose:
        Calculate fit to residual sky frequency. Sets attributes
        __f_sky_resid_fit, and updates __f_offset_fit with the new residual
        frequency fit

        Args:
            k (int): Order of polynomial fit made to residual frequency
            rho_exclude (list): Set of radius regions to exclude when making
                fit to residual frequency. Specify in km. Default is to
                exclude B ring region
            spm_include (list): Set of SPM regions to include when making fit
                to residual frequency. By default, only rho_exclude is used.
                If this keyword is specified, it overrides anything input to
                rho_exclude keyword. This is meant as an optional replacement
                of rho_exclude
            USE_GUI (bool): Use the interactive GUI to make a residual
                frequency fit. This is highly recommended
            TEST (bool): Optional test plot

        Dependencies:
            [1] RSRReader
            [2] Geometry
            [3] FResidFitGui
            [4] numpy
            [5] os
            [6] platform
            [7] time

        Notes:
        [1] Here is an example spm_include input, where each separate bracket
        contains a region that we want to include when fitting:
            spm_include = [[30250, 32600], [33520, 33890], [33990, 40200]]
        """

        # Attributes keeping track of input
        self._k = k
        self.__USE_GUI = USE_GUI

        # By default, exclude radii inside of Saturn, the B ring, and far
        #     outside of rings where sometimes the signal drops off. Only these
        #     regions are overly noisy at the default 8.192 second spacing in
        #     the frequency offset from FreqOffset class. This is the initial
        #     guess for regions to fit over. User adjusts using spm_include
        #     keyword
        #if rho_exclude is None:
        rho_exclude = [[0, 70000], [91900, 94000], [98000, 118000],
            [194400, np.inf]]

        f_rho = self.__f_rho
        f_spm = self.__f_spm

        # Residual frequency is amount of frequency offset not accounted
        # for by error in predicted spacecraft trajectory
        f_sky_resid = self.__f_offset - (self.__f_sky_recon - self.__f_sky_pred)

        # Determine if Cassini blocked by Saturn's ionosophere. Important for
        #     chord occultations
        is_blocked_atm, is_blocked_ion = cassini_blocked(f_spm, self.__rsr_inst,
            self.__kernels)

        # Array of indices to include by default. Overridden by spm_include
        #     keyword
        ind = []
        for i in range(len(rho_exclude)-1):
            ind.append(np.argwhere((f_rho > rho_exclude[i][1]) &
                (f_rho < rho_exclude[i+1][0]) &
                (np.invert(is_blocked_ion))))
        ind = np.reshape(np.concatenate(ind), -1)

        # If user specified spm_include argument that overrides the rho_exclude
        #     argument
        if spm_include is not None:
            self._spm_include = spm_include
            ind = []
            for i in range(len(spm_include)):
                ind.append(np.argwhere((f_spm >= spm_include[i][0]) &
                    (f_spm <= spm_include[i][1]) &
                    (np.invert(is_blocked_ion))))
            ind = np.reshape(np.concatenate(ind), -1)
        else:
            spm_include = []
            for i in range(len(rho_exclude)-1):
                _ind = np.argwhere((f_rho > rho_exclude[i][1])
                    & (f_rho < rho_exclude[i+1][0]))
                spm_include.append([float(min(f_spm[_ind])),
                    float(max(f_spm[_ind]))])
            self._spm_include = spm_include

        if len(ind) == 0:
            print('ERROR (FreqOffsetFit.set_f_sky_resid_fit): All specified \n'
                +'points fall withint specified rho exclusion zone')
            sys.exit()

        # When fitting, use x values adjusted to range over [-1, 1]
        npts = len(f_spm)
        spm_temp = (f_spm - f_spm[int(npts/2)])/max(f_spm - f_spm[int(npts/2)])

        # Coefficients for least squares fit, and evaluation of coefficients
        coef = poly.polyfit(spm_temp[ind], f_sky_resid[ind], k)
        f_sky_resid_fit = poly.polyval(spm_temp, coef)

        self.__f_sky_resid_fit = f_sky_resid_fit

        # Update frequency offset fit attribute for new residual sky
        # frequency fit
        self.__f_offset_fit = self.__f_sky_resid_fit + (self.__f_sky_recon -
            self.__f_sky_pred)

        if USE_GUI:
            root = Tk()
            f_resid_fit_gui_inst = FResidFitGui(root, self, f_spm, f_sky_resid)
            root.geometry('900x700+500+150')
            while True:
                try:
                    root.mainloop()
                    break
                except UnicodeDecodeError:
                    pass

        # Set history attribute
        self.__set_history()


    def get_f_sky_pred(self):
        """Returns predicted sky frequency and SPM it was evaluated at"""
        return self.__f_spm, self.__f_sky_pred


    def get_f_sky_recon(self):
        """Returns reconstructed sky frequency and SPM it was evaluated at"""
        return self.__f_spm, self.__f_sky_recon


    def get_f_sky_resid_fit(self):
        """Returns fit to residual frequency"""
        return self.__f_spm, self.__f_sky_resid_fit


    def get_f_offset_fit(self):
        """Returns fit to frequency offset"""
        return self.__f_spm, self.__f_offset_fit


    def get_IQ_c(self, spm_vals=None, IQ_m=None):
        """
        Purpose:
        Apply frequency offset fit to raw measured signal. Can supply
        SPM and IQ_m input if desired, otherwise the default is raw
        resolution

        Args:
            spm_vals (np.ndarray): SPM values corresponding to complex signal
                you're frequency correcting. Default is raw resolution from
                rsr_inst
            IQ_m (np.ndarray): Complex signal you're frequency correcting.
                Default is raw resolution from rsr_inst

        Outputs:
            spm_vals (np.ndarray): SPM values of the returned frequency
                corrected complex signal
            IQ_c (np.ndarray): Frequency corrected complex signal

        Dependencies:
            [1] numpy
            [2] scipy.interpolate

        Warnings:
            [1] If you resample IQ_m before frequency correcting it, you risk
                making the unnormalized signal much lower than it would
                otherwise be. Highly recommended to just use the default raw
                resolution spm_vals and IQ_m
            """

        if IQ_m is None:
            spm_vals = self.__spm_vals
            IQ_m = self.__IQ_m

        if len(spm_vals) != len(IQ_m):
            sys.exit('ERROR (FreqOffsetFit.get_IQ_c):'+
                     '\n SPM and IQ_m input lengths don\'t match')

        f_offset_fit = self.__f_offset_fit

        # Interpolate frequeny offset fit to 0.1 second spacing, since
        # this makes the integration later more accurate
        dt = 0.1
        npts = round((self.__f_spm[-1] - self.__f_spm[0])/dt)
        f_spm_interp = self.__f_spm[0] + dt*np.arange(npts)
        f_offset_fit_function = interp1d(self.__f_spm, f_offset_fit,
            fill_value='extrapolate')
        f_offset_fit_interp = f_offset_fit_function(f_spm_interp)

        # Integration of frequency offset fit to get phase detrending function.
        # Then interpolated to same SPM as I and Q
        f_detrend_interp = np.cumsum(f_offset_fit_interp)*dt
        f_detrend_interp_rad = f_detrend_interp*spice.twopi()
        f_detrend_rad_function = interp1d(f_spm_interp, f_detrend_interp_rad,
            fill_value='extrapolate')
        f_detrend_rad = f_detrend_rad_function(spm_vals)

        # Apply detrending function
        IQ_c = IQ_m*np.exp(-1j*f_detrend_rad)

        return spm_vals, IQ_c


    def __set_history(self):
        """
        Record information about the run and how to reproduce it.
        """

        # TODO (gsteranka): Replace entry to 'geo_inst' category of
        #     input_var_dict with self.__geo_inst.history
        input_var_dict = {'rsr_inst': self.__rsr_inst.history,
            'geo_inst': self.__geo_inst.history, 'f_uso': self.__f_uso}
        input_kw_dict = {'k': self._k, 'spm_include': self._spm_include,
            'USE_GUI': self.__USE_GUI}
        hist_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__,
            'Input Variables': input_var_dict,
            'Input Keywords':input_kw_dict}
        self.history = hist_dict


if __name__ == '__main__':
    pass
