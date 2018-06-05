#!/usr/bin/env python
"""

power_normalization_v2.py

Purpose: Normalize frequency-corrected power using a spline fit of specified
         order.

NOTE: Dependent on final format of geometry instances ("geo_inst")

Revisions:
      power_normalization.py
   2018 Mar 20 - gsteranka - Edited to get rid of debug steps
   2018 Mar 20 - gsteranka - Fixed bugs with spm_vals_down getting defined
                             wrong, and with default knots sometimes stretching
                             beyond the input times to make a spline fit for.
                             Added error check if no input points fall within
                             specified freespace regions
   2018 Apr 04 - gsteranka - In initial downsampling step before resample_poly,
                             downsample IQ_c_raw instead of p_obs_raw. This
                             accompanies the recent change in norm_diff_class.py
   2018 Apr 04 - gsteranka - In loop to make knots_spm and knots_rho, change i
                             index to int(i)
   2018 Apr 05 - gsteranka - Added index sorting right before splrep to take
                             ingress files, since splrep expects x values to be
                             in order, which isn't true for ingress
   2018 Apr 25 - gsteranka - Added testing plot in get_spline_fit() method.
                             Changed defining loop of knots_spm and knots_rho
                             to loop over range(n_knots_km) instead of
                             knots_km_where
      power_normalization_v2.py
   2018 May 03 - gsteranka - Copied from prev. version. Includes rsr_inst as
                             input, which is used to determine if Cassini is
                             blocked by Saturn or its ionosphere. Uses dsn,
                             year, and doy attributes
   2018 May 09 - gsteranka - Added freespace_spm and knots_spm keywords
   2018 May 11 - gsteranka - Added _freespace_spm and _knots_spm attributes
   2018 Jun 04 - gsteranka - Added history attribute and __set_history() method
"""

import numpy as np
import os
import platform
from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import splrep
from scipy.interpolate import splev
import sys
import time

from tkinter import Tk

try:
    from cassini_blocked import cassini_blocked
    from power_fit_gui import PowerFitGui
except ImportError:
    from ..tools.cassini_blocked import cassini_blocked
    from .power_fit_gui import PowerFitGui

class Normalization(object):
    """
    Class to get a power normalizing spline fit for frequency corrected
    data uniformly spaced in time at raw resolution.

    Example:
        >>> fit_inst = FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset, \
                f_uso, kernels, k=k, rho_exclude=rho_exclude)
        >>> (spm_raw, IQ_c_raw) = fit_inst.get_IQ_c()
        >>> norm_inst = Normalization(spm_raw, IQ_c_raw, geo_inst)
        >>> spm_fit = np.arange(spm_raw[0], spm_raw[-1], 1.0)
        >>> spline_fit = norm_inst.get_spline_fit(spm_fit, k=k, \
                knots_km=knots_km, dt_down=dt_down, \
                freespace_km=freespace_km, TEST=TEST)

    Attributes:
        _dt_down (float): Time spacing to downsample to before making
                a spline fit
        _freespace_km (list): Set of radius values, in km, to treat as
                free space. These are the regions that the spline fit is made
                to. Specify in km. Be sure to include a region for each knot
                specified
        _freespace_spm (list): SPM version of _freespace_km. Use this one if
            it's not None. Overrides _freespace_km
        __IQ_c_raw (np.ndarray): Raw resolution frequency corrected I and Q
        _k (int): Order of the spline fit. Default order is 2
        _knots_km (list): List of knots for the spline fit. Should only
                be at places where there is data. Specify in km
        _knots_spm (list): SPM versino of _knots_km. Use this one if it's not
            None. Overrides _knots_km
        __rho_interp_func (scipy.interpolate.interpolate.interp1d): interp1d
            function to get set of rho values in km for a set of SPM
        __rsr_inst: Instance of the RSRReader class
        kernels (list): List of kernels used for geometry
        __spm_raw (np.ndarray): Raw resolution SPM values
        _spm_fit (np.ndarray): SPM values that spline fit was evaluated at
        history (dict): Recorded information about the run
    """

    def __init__(self, spm_raw, IQ_c_raw, geo_inst, rsr_inst, TEST=False):
        """
        Purpose:
        Instantiation defines raw resolution SPM, frequency corrected I
        and Q, and a function to interpolate radius to any SPM value. These
        are set as attributes

        Args:
            spm_raw (np.ndarray): Raw resolution array of SPM values
            IQ_c_raw (np.ndarray): Frequency corrected complex signal at
                raw resolution
            geo_inst: Instance of geometry class. Contains attributes
                t_oet_spm and rho_km_vals. Can create mock version from a
                geometry file using geo_file_into_instance.py
            TEST (bool): Optional boolean argument that, if True, prints out
                intermediate values

        Dependencies:
            [1] RSRReader
            [2] Geometry
            [3] numpy
            [4] scipy
            [5] os
            [6] platform
            [7] sys
            [8] time

        Warnings:
            [1] If IQ signal is not properly frequency-corrected (i.e. if the
                residual frequency fit from the FreqOffsetFit class is bad),
                then you will get problems in this routine, since it
                downsamples the signal.
        """

        self.__rsr_inst = rsr_inst
        self.__geo_inst = geo_inst
        self.kernels = geo_inst.kernels

        self.__spm_raw = spm_raw
        self.__IQ_c_raw = IQ_c_raw

        # Function to interpolate spm to rho values in km
        rho_interp_func = interp1d(np.asarray(geo_inst.t_oet_spm_vals),
            np.asarray(geo_inst.rho_km_vals), fill_value='extrapolate')
        self.__rho_interp_func = rho_interp_func

        # Default arguments for spline fit. Chosen because:
        #     k: Order higher than 2 gives too much freedom to spline fit.
        #     knots_km: By experimentation, these tend to give the most
        #               consistently good spline fits. Sometimes need to
        #               adjust if a file's free space power cuts off sooner
        #               or begins later than expected
        #     freespace_km: These regions are always free space and
        #                   uninterrupted by diffraction pattern or an
        #                   eccentric ringlet. These accompany the default
        #                   knots
        self._k = 2
        #self._knots_km = np.array([70445., 87400., 117730., 119950., 133550.,
        #    194269.])
        self._knots_km = [70445., 87400., 117730., 119950., 133550., 194269.]
        self._dt_down = 0.5
        self._freespace_km = [[69100.0, 73500], [87350.0, 87450.0],
            [117690.0, 117780.0], [119850.0, 120020.0],
            [133500.0, 133650.0], [137000.0, 194400.0]]

        # Start off using knots and freespace in terms of radius, then use
        #     these only if they're specified
        self._knots_spm = None
        self._freespace_spm = None

        # Evaluate spline fit over geometry instance spm values by default
        self._spm_fit = geo_inst.t_oet_spm_vals

        if TEST:
            sample_rate_hz = int(round(1.0/(self.__spm_raw[1] -
                self.__spm_raw[0])))
            print('\nGeometry SPM, rho:')
            for i in range(10):
                print(geo_inst.t_oet_spm_vals[i], geo_inst.rho_km_vals[i])
            print('\nRaw resolution SPM, rho at 1 sec. interval:')
            for i in range(10):
                print(self.__spm_raw[i*sample_rate_hz],
                    self.__rho_interp_func(self.__spm_raw[i*sample_rate_hz]))

        self.__set_history()


    def __downsample_IQ(self, dt_down, TEST=False):
        """
        Downsample complex signal to specified time spacing to avoid
        diffraction pattern ruining spline fit

        Args:
            dt_down (float): Time spacing to downsample to
            TEST (bool): If True, prints downsampled results

        Returns:
            spm_vals_down (np.ndarray): SPM values after downsampling
            rho_km_vals_down (np.ndarray): Rho values after downsampling
            p_obs_down (np.ndarray): Observed power after downsampling"""

        # Downsampling coefficient q
        dt_raw = self.__spm_raw[1] - self.__spm_raw[0]
        q = round(dt_down / dt_raw)

        # Downsample IQ_c by factor of q and not power because this is to
        #     match the resampling done in norm_diff_class.py, where it also
        #     resamples IQ_c
        IQ_c_down = signal.resample_poly(self.__IQ_c_raw, 1, q)
        p_obs_down = abs(IQ_c_down)**2

        # New SPM, rho, and p_obs, at the downsampled resolution
        spm_vals_down = np.linspace(self.__spm_raw[0], self.__spm_raw[-1],
            num=len(p_obs_down))
        rho_km_vals_down = self.__rho_interp_func(spm_vals_down)

        if TEST:
            print('\nDownsampled SPM, rho:')
            for i in range(10):
                print(spm_vals_down[i], rho_km_vals_down[i])

        return spm_vals_down, rho_km_vals_down, p_obs_down


    def get_spline_fit(self, spm_fit=None, k=None, knots_km=None,
            dt_down=None, freespace_km=None, freespace_spm=None,
            knots_spm=None, USE_GUI=True, TEST=False):
        """
        Purpose:
        Make spline fit to observed downsampled power at specified set
        of times spm_fit. Specified keyword arguments will override the
        corresponding defaults set in __init__

        Args:
            spm_fit (np.ndarray): SPM values to evaluate the spline fit at
            k (int): Order of the spline fit. Default order is 2
            knots_km (list): List of knots for the spline fit. Should only
                be at places where there is data. Specify in km
            dt_down (float): Time spacing to downsample to before making
                a spline fit
            freespace_km (list): Set of radius values, in km, to treat as
                free space. These are the regions that the spline fit is made
                to. Specify in km. Be sure to include a region for each knot
                specified
            freespace_spm (list): Set of SPM values to treat as free space.
                Meant as an optional replacement for freespace_km. Setting this
                will override the freespace_km keyword
            knots_spm (list): List of knots for the spline fit. Same as
                knots_km, but in SPM instead of radius. Specifying this
                overrides knots_km
            USE_GUI (bool): Use the interactive GUI to make a spline fit to
                power. This is highly recommended
            TEST (bool): If True, print out intermediate values

        Outputs:
            spm_fit: SPM values for the spline_fit output
            spline_fit: Spline fit to observed power at the specified times
                spm_fit

        Dependencies:
            [1] PowerFitGui
            [2] cassini_blocked
            [3] numpy
            [4] scipy
            [5] sys
            [6] os
            [7] platform
            [8] time
        """

        # Update defaults if any keyword arguments were specified
        if spm_fit is not None:
            self._spm_fit = spm_fit
        if k is not None:
            self._k = k
        if knots_km is not None:
            self._knots_km = knots_km
        if dt_down is not None:
            self._dt_down = dt_down
        if freespace_km is not None:
            self._freespace_km = freespace_km
        if freespace_spm is not None:
            self._freespace_spm = freespace_spm
        if knots_spm is not None:
            self._knots_spm = knots_spm

        # Use whatever most up-to-date arguments are
        spm_fit = self._spm_fit
        k = self._k
        knots_km = self._knots_km
        dt_down = self._dt_down
        freespace_km = self._freespace_km
        freespace_spm = self._freespace_spm
        knots_spm = self._knots_spm

        # Downsample I and Q to the time spacing dt_down
        (spm_vals_down, rho_km_vals_down,
         p_obs_down) = self.__downsample_IQ(dt_down, TEST=TEST)

        # Determine if Cassini is behind Saturn
        is_blocked_atm, is_blocked_ion = cassini_blocked(spm_vals_down,
            self.__rsr_inst, self.kernels)

         # Define freespace regions
        ind_free = []
        for i in range(len(freespace_km)):
            ind_free.append(np.argwhere((rho_km_vals_down > freespace_km[i][0])
                & (rho_km_vals_down < freespace_km[i][1])
                & (np.invert(is_blocked_atm))))
        ind_free = np.reshape(np.concatenate(ind_free), -1)

        if freespace_spm is not None:
            ind_free = []
            for i in range(len(freespace_spm)):
                ind_free.append(np.argwhere((spm_vals_down > freespace_spm[i][0])
                    & (spm_vals_down < freespace_spm[i][1])
                    & (np.invert(is_blocked_atm))))
            ind_free = np.reshape(np.concatenate(ind_free), -1)

        # Limit downsampled data to freespace regions
        spm_vals_free = spm_vals_down[ind_free]
        rho_km_vals_free = rho_km_vals_down[ind_free]
        p_obs_free = p_obs_down[ind_free]

        # Find knots in range of specified radii
        if knots_spm is None:
            knots_km_where = np.argwhere(
                np.logical_and(knots_km>min(rho_km_vals_free),
                    knots_km<max(rho_km_vals_free)))
            
            # Can't make fit if no specified knots lie in freespace regions
            if len(knots_km_where) == 0:
                print('ERROR (Normalization.get_spline_fit): No specified '
                    + 'knots \n'
                    + 'fall within freespace range. Freespace ranges are '
                    + str(freespace_km))
                sys.exit()

            # Select data points closest to selected knots
            n_knots_km = len(knots_km_where)
            knots_spm_data = np.zeros(n_knots_km)
            knots_rho_data = np.zeros(n_knots_km)
            for i in range(n_knots_km):
                _knot_km = knots_km[int(knots_km_where[i])]
                _ind = np.argmin(abs(rho_km_vals_free - _knot_km))
                knots_spm_data[i] = spm_vals_free[_ind]
                knots_rho_data[i] = rho_km_vals_free[_ind]

        # Find knots in range of specified SPM
        else:
            knots_spm_where = np.argwhere(
                np.logical_and(knots_spm>min(spm_vals_free),
                    knots_spm<max(spm_vals_free)))

            # Can't make fit if no specified knots lie in freespace regions
            if len(knots_spm_where) == 0:
                print('ERROR (Normalization.get_spline_fit): No specified '
                    + 'times \n'
                    + 'fall within freespace range. reespace ranges are '
                    + str(freespace_spm))
                sys.exit()

            # Select data points closest to selected knots
            n_knots_spm = len(knots_spm_where)
            knots_spm_data = np.zeros(n_knots_spm)
            knots_rho_data = np.zeros(n_knots_spm)
            for i in range(n_knots_spm):
                _knot_spm = knots_spm[int(knots_spm_where[i])]
                _ind = np.argmin(abs(spm_vals_free - _knot_spm))
                knots_spm_data[i] = spm_vals_free[_ind]
                knots_rho_data[i] = rho_km_vals_free[_ind]

        if TEST:
            print('\nSelected knots:')
            print(knots_spm_data)
            print(knots_rho_data)

        # Make and evaluate the spline fit. Sort spm_vals_free and p_obs_free
        #     in case it's ingress, since indices were gotten from rho values,
        #     which are not in order for ingress
        ind_sort = np.argsort(spm_vals_free)
        ind_knot_sort = np.argsort(knots_spm_data)
        spline_rep = splrep(spm_vals_free[ind_sort], p_obs_free[ind_sort],
            k=k, t=knots_spm_data[ind_knot_sort])
        spline_fit = splev(spm_fit, spline_rep)

        if USE_GUI:
            root = Tk()
            power_fit_gui_inst = PowerFitGui(root, self, spm_vals_down,
                p_obs_down, spm_fit)
            root.geometry('900x700+500+150')
            while True:
                try:
                    root.mainloop()
                    break
                except UnicodeDecodeError:
                    pass
            spline_fit = power_fit_gui_inst.yfit
            self._k = power_fit_gui_inst.fit_deg
            self._knots_spm = power_fit_gui_inst.knots_spm
            self._freespace_spm = power_fit_gui_inst.xlim

        self.__set_history()
        return spm_fit, spline_fit


    def __set_history(self):
        """
        Record info about the run's history
        """

        input_var_dict = {'spm_raw': self.__spm_raw,
            'IQ_c_raw': self.__IQ_c_raw, 'geo_inst': self.__geo_inst.history,
            'rsr_inst': self.__rsr_inst.history}
        input_kw_dict = {'spm_fit': self._spm_fit, 'k': self._k,
            'knots_km': self._knots_km, 'dt_down': self._dt_down,
            'freespace_km': self._freespace_km,
            'freespace_spm': self._freespace_spm, 'knots_spm': self._knots_spm}
        hist_dict = {'user name': os.getlogin(),
            'host name': os.uname().nodename,
            'run date': time.ctime() + ' ' + time.tzname[0],
            'python version': platform.python_version(),
            'operating system': os.uname().sysname,
            'source file': __file__,
            'input variables': input_var_dict,
            'input keywords':input_kw_dict}
        self.history = hist_dict


if __name__ == '__main__':
    pass
