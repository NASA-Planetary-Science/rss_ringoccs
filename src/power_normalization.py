#!/usr/bin/env python
"""

power_normalization.py

Purpose: Normalize frequency-corrected power using a spline fit of specified
         order.

NOTE: Dependent on final format of geometry instances ("geo_inst")

Revisions:
      gjs_power_normalization.py
   Feb 26 2018 - gsteranka - Original version
      gjs_power_normalization_v2.py
   Mar 06 2018 - gsteranka - Edited so user only needs to instantiate and
                             run get_spline_fit, which now returns spline fit
                             for a specified array of SPM values.
      gjs_power_normalization_v3.py
   Mar 19 2018 - gsteranka - Edited so default inputs are set as attributes in
                             __init__, and the inputs in get_spline_fit can
                             override those defaults as the new attribute
                             values. Done this way to keep track of what the
                             previous fit specifications were inside of the
                             instance.
      power_normalization.py
   Mar 20 2018 - gsteranka - Edited to get rid of debug steps
   Mar 20 2018 - gsteranka - Fixed bugs with spm_vals_down getting defined
                             wrong, and with default knots sometimes stretching
                             beyond the input times to make a spline fit for.
                             Added error check if no input points fall within
                             specified freespace regions
   Apr 04 2018 - gsteranka - Downsample IQ_c_raw instead of p_obs_raw. This
                             accompanies the recent change in
                             norm_diff_class.py
   Apr 04 2018 - gsteranka - In loop to make knots_spm and knots_rho, change i
                             index to int(i)
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import splrep
from scipy.interpolate import splev
import sys

class Normalization(object):
    """Class to get a power normalizing spline fit for frequency corrected
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
        _k (int): Order of the spline fit. Default order is 2
        _knots_km (list): List of knots for the spline fit. Should only
                be at places where there is data. Specify in km
        _dt_down (float): Time spacing to downsample to before making
                a spline fit
        _freespace_km (list): Set of radius values, in km, to treat as
                free space. These are the regions that the spline fit is made
                to. Specify in km. Be sure to include a region for each knot
                specified
        __spm_raw (np.ndarray): Raw resolution SPM values
        __IQ_c_raw (np.ndarray): Raw resolution frequency corrected I and Q
        __rho_interp_func (scipy.interpolate.interpolate.interp1d): interp1d
            function to get set of rho values in km for a set of SPM
    """

    def __init__(self, spm_raw, IQ_c_raw, geo_inst, TEST=False):
        """Instantiation defines raw resolution SPM, frequency corrected I
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
                intermediate values"""

        self.__spm_raw = spm_raw
        self.__IQ_c_raw = IQ_c_raw

        # Function to interpolate spm to rho values in km
        rho_interp_func = interp1d(np.asarray(geo_inst.t_oet_spm_vals),
                                   np.asarray(geo_inst.rho_km_vals),
                                   fill_value='extrapolate')
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
        self._knots_km = np.array([70445., 87400., 117730., 119950., 133550.,
                                   194269.])
        self._dt_down = 0.5
        self._freespace_km = [[69100.0, 73500], [87350.0, 87450.0],
                             [117690.0, 117780.0], [119850.0, 120020.0],
                             [133500.0, 133650.0], [137000.0, 194400.0]]

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


    def __downsample_IQ(self, dt_down, TEST=False):
        """Downsample complex signal to specified time spacing to avoid
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
        # match the resampling done in norm_diff_class.py, where it also
        # resamples IQ_c
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


    def get_spline_fit(self, spm_fit, k=None, knots_km=None,
                       dt_down=None, freespace_km=None, TEST=False):
        """Make spline fit to observed downsampled power at specified set
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
            TEST (bool): If True, print out intermediate values

        Returns:
            spline_fit: Spline fit to observed power at the specified times
                spm_fit"""

        # Update defaults if any keyword arguments were specified
        if k is not None:
            self._k = k
        if knots_km is not None:
            self._knots_km = knots_km
        if dt_down is not None:
            self._dt_down = dt_down
        if freespace_km is not None:
            self._freespace_km = freespace_km

        # Use whatever most up-to-date arguments are
        k = self._k
        knots_km = self._knots_km
        dt_down = self._dt_down
        freespace_km = self._freespace_km

        # Downsample I and Q to the time spacing dt_down
        (spm_vals_down, rho_km_vals_down,
         p_obs_down) = self.__downsample_IQ(dt_down, TEST=TEST)

         # Define freespace regions
        ind_free = []
        for i in range(len(freespace_km)):
            ind_free.append(np.argwhere((rho_km_vals_down > freespace_km[i][0])
                                 & (rho_km_vals_down < freespace_km[i][1])))
        ind_free = np.reshape(np.concatenate(ind_free), -1)

        # Limit downsampled data to freespace regions
        spm_vals_free = spm_vals_down[ind_free]
        rho_km_vals_free = rho_km_vals_down[ind_free]
        p_obs_free = p_obs_down[ind_free]

        # Find knots in range of specified times
        knots_km_where = np.argwhere(
            np.logical_and(knots_km>min(rho_km_vals_free),
                           knots_km<max(rho_km_vals_free)))

        # If specified times don't have freespace region, you can't make a fit
        if len(knots_km_where) == 0:
            print('ERROR (Normalization.get_spline_fit): No specified times \n'+
                  'fall within freespace range. Freespace ranges are '+
                  str(freespace_km))
            sys.exit()

        # Select data points closest to selected knots
        n_knots_km = len(knots_km_where)
        knots_spm = np.zeros(n_knots_km)
        knots_rho = np.zeros(n_knots_km)
        for i in knots_km_where:#range(n_knots_km):
            _ind = np.argmin(abs(rho_km_vals_free - knots_km[int(i)]))
            knots_spm[i] = spm_vals_free[_ind]
            knots_rho[i] = rho_km_vals_free[_ind]

        if TEST:
            print('\nSelected knots:')
            print(knots_spm)
            print(knots_rho)

        # Make and evaluate the spline fit. Sort spm_vals_free and p_obs_free
        # in case it's ingress, since indices were gotten from rho values,
        # which are not in order for ingress
        ind_sort = np.argsort(spm_vals_free)
        ind_knot_sort = np.argsort(knots_spm)
        spline_rep = splrep(spm_vals_free[ind_sort], p_obs_free[ind_sort],
                            k=k, t=knots_spm[ind_knot_sort])
        spline_fit = splev(spm_fit, spline_rep)

        return spline_fit


if __name__ == '__main__':
    pass
