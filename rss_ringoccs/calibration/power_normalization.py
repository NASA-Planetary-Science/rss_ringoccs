#!/usr/bin/env python
"""

power_normalization.py

Purpose: Normalize frequency-corrected power using a spline fit of specified
         order.

WHERE TO GET NECESSARY INPUT:
    spm_raw, IQ_c_raw: Use the output of the "get_IQ_c" method in the
        FreqOffsetFit class, found in
        rss_ringoccs/calibration/freq_offset_fit.py
    geo_inst: Use an instance of the Geometry class, found inside of
        rss_ringoccs/occgeo/calc_occ_geometry.py
    rsr_inst: Use an instance of the RSRReader class, found inside of
        rss_ringoccs/rsr_reader/rsr_reader.py

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
                             accompanies the recent change in
                             norm_diff_class.py
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
   2018 Sep 19 - jfong - removed print statement of spm and rho values,
                            and for selected knots
    2018 Sep 20 - jfong - add file_search kwarg, search file for splrep outputs
                        - add profdir attribute
    2018 Sep 25 - jfong - remove file_search kwarg from __init__ because
                          get_spline_fit() is never called within __init__ but
                          always outside of the function (which needs to be
                          addressed later)
                        - removed try/except for faulty spline orders, add
                          the spline order check when extracting order from
                          kwarg
    2018 Sep 26 - jfong - write PNFP file if file not found,
                          read PNFP file if specified in file_search
"""

import numpy as np
import pdb
from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import splrep
from scipy.interpolate import splev
import sys
import pickle

from tkinter import Tk

sys.path.append('../..')
import rss_ringoccs as rss
from rss_ringoccs.tools.search_for_file import search_for_file
from rss_ringoccs.tools.write_intermediate_files import write_intermediate_files
from rss_ringoccs.tools.write_output_files import write_output_files
sys.path.remove('../..')

from ..tools.cassini_blocked import cassini_blocked
from .power_fit_gui import PowerFitGui



class Normalization(object):
    """
    Class to get a power normalizing spline fit for frequency corrected
    data uniformly spaced in time at raw resolution.

    Example:
        >>> rsr_inst = rss.rsr_reader.RSRReader(rsr_file)
        >>> geo_inst = rss.occgeo.Geometry(rsr_inst, 'Saturn', 'Cassini',
                kernels)
        >>> fit_inst = rss.calibration.FreqOffsetFit(rsr_inst, geo_inst, f_spm,
                f_offset, f_uso, kernels, poly_order=poly_order,
                spm_include=spm_include)
        >>> spm_raw, IQ_c_raw = fit_inst.get_IQ_c()
        >>> norm_inst = rss.calibration.Normalization(spm_raw, IQ_c_raw,
                geo_inst, rsr_inst)
        >>> spm_fit, spline_fit = norm_inst.get_spline_fit(
                spline_order=spline_order, knots_spm=knots_spm,
                dt_down=dt_down, freespace_spm=freespace_spm, verbose=verbose)

    Attributes:
        _dt_down (float):
            Time spacing to downsample to before making a spline fit
        _freespace_km (list):
            Set of default radius values, in km, to treat as
            free space. These are the regions that the spline fit is made
            to. Specify in km. Be sure to include a region for each knot
            specified
        _freespace_spm (list):
            SPM version of _freespace_km. Use this one if
            it's not None. Overrides _freespace_km
        __IQ_c_raw (np.ndarray):
            Raw resolution frequency corrected I and Q
        _spline_order (int):
            Order of the spline fit. Default order is 2
        _knots_km (list):
            List of default knot positions for the spline fit
        _knots_spm (list):
            SPM versino of _knots_km. Use this one if it's not
            None. Overrides _knots_km
        __rho_interp_func (scipy.interpolate.interpolate.interp1d):
            interp1d function to get set of rho values in km for a set of SPM
        __rsr_inst:
            Instance of the RSRReader class
        kernels (list):
            List of kernels used for geometry
        __spm_raw (np.ndarray):
            Raw resolution SPM values
        _spm_fit (np.ndarray):
            SPM values that spline fit was evaluated at
        history (dict):
            Recorded information about the run
    """

    def __init__(self, spm_raw, IQ_c_raw, geo_inst, rsr_inst, verbose=False):
        """
        Purpose:
        Instantiation defines raw resolution SPM, frequency corrected I
        and Q, and a function to interpolate radius to any SPM value. These
        are set as attributes

        Args:
            spm_raw (np.ndarray):
                Raw resolution array of SPM values. Get from
                using get_IQ_c() method of the FreqOffsetFit class
            IQ_c_raw (np.ndarray):
                Frequency corrected complex signal at
                raw resolution. Get from using get_IQ_c() method of the
                FreqOffsetFit class
            geo_inst:
                Instance of geometry class. Contains attributes
                t_oet_spm and rho_km_vals. Can create mock version from a
                geometry file using geo_file_into_instance.py
            verbose (bool):
                Optional boolean argument that, if True, prints out
                intermediate values

        Dependencies:
            [1] RSRReader
            [2] Geometry
            [3] numpy
            [4] scipy
            [5] sys

        Warnings:
            [1] If IQ signal is not properly frequency-corrected (i.e. if the
                residual frequency fit from the FreqOffsetFit class is bad),
                then you will get problems in this routine, since it
                downsamples the signal.
        """

        if len(spm_raw) != len(IQ_c_raw):
            print('ERROR (Normalization): Input arrays not the same length. '
                + 'Input should be from the get_IQ_c() method of the '
                + 'FreqOffsetFit class')
        elif (type(spm_raw) != np.ndarray) | (type(IQ_c_raw) != np.ndarray):
            print('ERROR (Normalization): Input arrays must be numpy arrays. '
                + 'Input should be from the get_IQ_c() method of the '
                + 'FreqOffsetFit class')

        if not isinstance(geo_inst, rss.occgeo.Geometry):
            sys.exit('ERROR (Normalization): geo_inst input must be an '
                + 'instance of the Geometry class')

        if not isinstance(rsr_inst, rss.rsr_reader.RSRReader):
            sys.exit('ERROR (Normalization): rsr_inst input must be an '
                + 'instance of the RSRReader class')

        if not isinstance(verbose, bool):
            print('WARNING (Normalization): verbose input must be boolean. '
                + 'Ignoring current input and setting to False')
            verbose = False

        self.__rsr_inst = rsr_inst
        self.__geo_inst = geo_inst
        self.profdir = geo_inst.get_profile_dir()
        self.kernels = geo_inst.kernels

        self.__spm_raw = spm_raw
        self.__IQ_c_raw = IQ_c_raw

        # Function to interpolate spm to rho values in km
        rho_interp_func = interp1d(np.asarray(geo_inst.t_oet_spm_vals),
            np.asarray(geo_inst.rho_km_vals), fill_value='extrapolate')
        self.__rho_interp_func = rho_interp_func

        # Default arguments for spline fit. Mentioned in docstring for
        #     get_spline_fit method, under the descriptions for the inputs
        #     spline_order, freespace_spm, and knots_spm
        self._spline_order = 2
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

        self.__set_history()

    def __downsample_IQ(self, dt_down, verbose=False):
        """
        Purpose:
        Downsample complex signal to specified time spacing to avoid
        diffraction pattern ruining spline fit

        Args:
            dt_down (float):
                Time spacing to downsample to
            verbose (bool):
                If True, prints downsampled results

        Returns:
            spm_vals_down (np.ndarray):
                SPM values after downsampling
            rho_km_vals_down (np.ndarray):
                Rho values after downsampling
            p_obs_down (np.ndarray):
                Observed power after downsampling"""

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

        return spm_vals_down, rho_km_vals_down, p_obs_down

    def get_spline_fit(self, spline_order=None, dt_down=None,
            freespace_spm=None, knots_spm=None, USE_GUI=True, verbose=False,
            file_search=False):
        """
        Purpose:
        Make spline fit to observed downsampled power at specified set
        of times self._spm_fit, which are equal to geometry times. Specified
        keyword arguments will override the corresponding defaults set in
        __init__

        Args:
            spline_order (int):
                Order of the spline fit. Default order is 2
                because any order higher than this gives the splien fit too
                much freedom
            dt_down (float):
                Time spacing to downsample to before making a spline fit
            freespace_spm (list):
                Set of SPM values to treat as free space.
                Meant as an optional replacement for the defaults specified by
                the _freespace_km attribute. Setting this will override the
                default. Defaults were chosen based on what works for rev7E
            knots_spm (list):
                List of knots for the spline fit in SPM (on the
                bottom of the GUI x axis). Specifying this overrides knots_km.
                The _knots_km  attribute gives the default knots, which were
                chosen based on what worked for rev7E. Specifying knots_spm
                keyword overrides the default.
            USE_GUI (bool):
                Use the interactive GUI to make a spline fit to
                power. This is highly recommended
            verbose (bool):
                If True, print out intermediate values
            file_search (bool, str):
                Search ../output/ directory for power normalization fit
                parameter (PNFP) file if set. If string given, read fit from
                that particular file. If not set, compute freespace power
                spline fit and save to PNFP file.

        Outputs:
            spm_fit:
                SPM values for the spline_fit output
            spline_fit:
                Spline fit to observed power at the specified times spm_fit

        Dependencies:
            [1] PowerFitGui
            [2] cassini_blocked
            [3] numpy
            [4] scipy
            [5] sys

        Notes:
            [1] HIGHLY RECOMMENDED to use the GUI if you haven't done the
                occultation before. Much easier to tinker with the fit this way

        Warnings:
            [1] If you make dt_down outrageously large, like 1000s, than
                you'll have a hard time making a good fit. If you make it
                longer than the length of the occultation, like 20000s,
                than you'll probably get an error
            [2] Setting file_search=False and USE_GUI=False will not save any
                intermediate files.
        """

        if not isinstance(USE_GUI, bool):
            print('WARNING (Normalization): USE_GUI input must be boolean. '
                  + 'Ignoring current input and setting to True')
            USE_GUI = True

        if not isinstance(verbose, bool):
            print('WARNING (Normalization): verbose input must be boolean. '
                  + 'Ignoring current input and setting to False')
            verbose = False

        # SPM values to evaluate spline fit at
        spm_fit = self._spm_fit

        # Update defaults if any keyword arguments were specified
        if spline_order is not None:
            if spline_order<1 or spline_order>5:
                print('WARNING (Normalization.get_spline_fit): Given degree of '
                        + 'the spline (k) is not supported (1<=k<=5). '
                        + 'Reverting to degree 2')
            else:
                self._spline_order = spline_order
        if dt_down is not None:
            dt_down = abs(dt_down)
            self._dt_down = dt_down
        if freespace_spm is not None:
            self._freespace_spm = freespace_spm
        if knots_spm is not None:
            self._knots_spm = knots_spm

        # Use whatever most up-to-date arguments are
        spline_order = self._spline_order
        knots_km = self._knots_km
        dt_down = self._dt_down
        freespace_km = self._freespace_km
        freespace_spm = self._freespace_spm
        knots_spm = self._knots_spm

        # Downsample I and Q to the time spacing dt_down
        if verbose:
            print('\tDownsampling input IQ_c to ' + str(dt_down)
                + 's time spacing...')
        try:
            (spm_vals_down, rho_km_vals_down,
                p_obs_down) = self.__downsample_IQ(dt_down, verbose=verbose)
        except TypeError:
            print('WARNING (Normalization.get_spline_fit): dt_down must be '
                + 'either int or float. Reverting to the default of 0.5s')
            self._dt_down = 0.5
            dt_down = 0.5
            (spm_vals_down, rho_km_vals_down,
                p_obs_down) = self.__downsample_IQ(dt_down, verbose=verbose)

        # Determine if Cassini is behind Saturn
        if verbose:
            print('\tFinding where atmosphere and ionosphere occults '
                + 'spacecraft...')
        is_blocked_atm, is_blocked_ion = cassini_blocked(spm_vals_down,
            self.__rsr_inst, self.kernels)

        # Define freespace regions. Overridden if freespace_spm is not None
        ind_free = self.__get_default_freespace_ind(rho_km_vals_down,
            freespace_km, is_blocked_atm)

        # Override default freespace values above
        try:
            if freespace_spm is not None:
                ind_free = []
                for i in range(len(freespace_spm)):
                    ind_free.append(np.argwhere(
                        (spm_vals_down > freespace_spm[i][0])
                        & (spm_vals_down < freespace_spm[i][1])
                        & (np.invert(is_blocked_atm))))
                ind_free = np.reshape(np.concatenate(ind_free), -1)
        except TypeError:
            print('WARNING (Normalization.get_spline_fit): Illegal input for '
                + 'freespace_spm. Reverting to default freespace regions')
            self._freespace_spm = None
            freespace_spm = None
            ind_free = self.__get_default_freespace_ind(rho_km_vals_down,
                freespace_km, is_blocked_atm)

        if np.all(ind_free is False):
            print('WARNING (Normalization.get_spline_fit): No SPM values fall '
                + 'within input freespace_spm regions. Returning to default '
                + 'freespace regions')
            self._freespace_spm = None
            freespace_spm = None
            ind_free = self.__get_default_freespace_ind(rho_km_vals_down,
                freespace_km, is_blocked_atm)

        # Limit downsampled data to freespace regions
        spm_vals_free = spm_vals_down[ind_free]
        rho_km_vals_free = rho_km_vals_down[ind_free]
        p_obs_free = p_obs_down[ind_free]

        # Find knots in range of specified radii
        if knots_spm is None:
            knots_km_where = np.argwhere(
                np.logical_and(knots_km > min(rho_km_vals_free),
                    knots_km < max(rho_km_vals_free)))

            # Can't make fit if no specified knots lie in freespace regions
            if len(knots_km_where) == 0:
                print('WARNING (Normalization.get_spline_fit): No default '
                    + 'knots fall within freespace range. Adjust freespace_spm'
                    + 'input. For now, returning to default freespace range')
                self._freespace_spm = None
                freespace_spm = None
                ind_free = self.__get_default_freespace_ind(rho_km_vals_down,
                    freespace_km, is_blocked_atm)
                spm_vals_free = spm_vals_down[ind_free]
                rho_km_vals_free = rho_km_vals_down[ind_free]
                p_obs_free = p_obs_down[ind_free]
                knots_km_where = np.argwhere(
                    np.logical_and(knots_km > min(rho_km_vals_free),
                        knots_km < max(rho_km_vals_free)))

            # Select data points closest to selected knots
            knots_spm_data, knots_rho_data = self.__get_ind_from_knots(
                'RHO_KM', knots_km_where, knots_km, spm_vals_free,
                rho_km_vals_free)

        # Find knots in range of specified SPM
        else:
            knots_spm_where = np.argwhere(
                np.logical_and(knots_spm > min(spm_vals_free),
                    knots_spm < max(spm_vals_free)))

            # Input for __get_ind_from_knots
            type_of_knot = 'SPM'
            knots_where_input = knots_spm_where
            knots_input = knots_spm

            # Can't make fit if no specified knots lie in freespace regions
            if len(knots_spm_where) == 0:
                print('WARNING (Normalization.get_spline_fit): No specified '
                    + 'knots fall within freespace range. Returning to '
                    + 'default knots. Freespace ranges are currently: '
                    + '\n\tfreespace_spm = ' + str(freespace_spm)
                    + '\n\tfreespace_km = ' + str(freespace_km))
                self._knots_spm = None
                knots_spm = None
                knots_km_where = np.argwhere(
                    np.logical_and(knots_km > min(rho_km_vals_free),
                        knots_km < max(rho_km_vals_free)))

                # Input for __get_ind_from_knots if no specified knots lie
                #     in freespace regions
                type_of_knot = 'RHO_KM'
                knots_where_input = knots_km_where
                knots_input = knots_km

            # Select data points closest to selected knots
            knots_spm_data, knots_rho_data = self.__get_ind_from_knots(
                type_of_knot, knots_where_input, knots_input, spm_vals_free,
                rho_km_vals_free)

        # Make and evaluate the spline fit. Sort spm_vals_free and p_obs_free
        #     in case it's ingress, since indices were gotten from rho values,
        #     which are not in order for ingress
        if verbose:
            print('\tEvaluating spline fit...')
        ind_sort = np.argsort(spm_vals_free)
        ind_knot_sort = np.argsort(knots_spm_data)
        if file_search:
            if isinstance(file_search, str):
                print('\tExtracting power normalization spline fit parameters '
                        + 'from:\n\t\t' + '/'.join(file_search.split('/')[0:5])
                        + '/\n\t\t\t' + file_search.split('/')[-1])
                file_object = open(file_search, 'rb')
                fit_param_dict = pickle.load(file_object)
                k_power_norm = fit_param_dict['spline_order']
                spline_coef = fit_param_dict['spline_coef']
                knots_spm = fit_param_dict['knots_spm']
                spline_rep = (knots_spm, spline_coef, k_power_norm)
            else:
                # Search for power normalization fit pickle (PNFP) file
                pnfp_file = search_for_file(self.__rsr_inst.year,
                        self.__rsr_inst.doy, self.__rsr_inst.band,
                        self.__rsr_inst.dsn, self.profdir, 'PNFP')
                if pnfp_file == 'N/A':
                    print('WARNING (get_spline_fit()): PNFP file not found!')
                    print('Evaluating spline fit...')
                    spline_rep = splrep(spm_vals_free[ind_sort],
                            p_obs_free[ind_sort],
                            k=spline_order, t=knots_spm_data[ind_knot_sort])
                    pnfp = {'knots_spm': spline_rep[0],
                            'spline_coef': spline_rep[1],
                            'spline_order': spline_rep[2]}
                    write_intermediate_files(self.__rsr_inst.year,
                            self.__rsr_inst.doy, self.__rsr_inst.band,
                            self.__rsr_inst.dsn, self.profdir, 'PNFP', pnfp)
                else:
                    print('\tExtracting power normalization spline fit parameters '
                            + 'from:\n\t\t' + '/'.join(pnfp_file.split('/')[0:5])
                            + '/\n\t\t\t' + pnfp_file.split('/')[-1])
                    file_object = open(pnfp_file, 'rb')
                    fit_param_dict = pickle.load(file_object)
                    k_power_norm = fit_param_dict['spline_order']
                    spline_coef = fit_param_dict['spline_coef']
                    knots_spm = fit_param_dict['knots_spm']
                    spline_rep = (knots_spm, spline_coef, k_power_norm)
        else:
            print('\tEvaluating spline fit...')
            spline_rep = splrep(spm_vals_free[ind_sort],
                    p_obs_free[ind_sort],
                    k=spline_order, t=knots_spm_data[ind_knot_sort])

        spline_fit = splev(spm_fit, spline_rep, ext=1)

        # Set values outside of SPM values in splrep to the exterior-most
        #     spline values
        try:
            min_fit_ind = np.max(((spline_fit == 0)
                & (spm_fit <= min(spm_vals_free[ind_sort]))).nonzero())
            spline_fit[0:min_fit_ind + 1] = (
                np.zeros(len(spline_fit[0:min_fit_ind + 1]))
                + spline_fit[min_fit_ind + 1])
        except ValueError:
            pass
        try:
            max_fit_ind = np.min(((spline_fit == 0)
                & (spm_fit >= max(spm_vals_free[ind_sort]))).nonzero())
            spline_fit[max_fit_ind:] = (np.zeros(len(spline_fit[max_fit_ind:]))
                + spline_fit[max_fit_ind - 1])
        except ValueError:
            pass

        self.pnfp_splrep = spline_rep

        if USE_GUI:
            print('Using GUI')
            pnfp_coef_ori = spline_rep[1]
            root = Tk()
            power_fit_gui_inst = PowerFitGui(root, self, spm_vals_down,
                rho_km_vals_down, p_obs_down, spm_fit)
            root.geometry('900x700+500+150')
            while True:
                try:
                    root.mainloop()
                    break
                except UnicodeDecodeError:
                    pass
            spline_fit = power_fit_gui_inst.yfit

            # If fit was changed during GUI, write new file
            if len(pnfp_coef_ori) == len(self.pnfp_splrep[1]):
                if all(pnfp_coef_ori == self.pnfp_splrep[1]) == False:
                    print('\tPower normalization fit changed within GUI!\n'
                            + '\tWriting new file...')
                    pnfp = {'knots_spm': self.pnfp_splrep[0],
                            'spline_coef': self.pnfp_splrep[1],
                            'spline_order': self.pnfp_splrep[2]}
                    write_intermediate_files(self.__rsr_inst.year,
                            self.__rsr_inst.doy, self.__rsr_inst.band,
                            self.__rsr_inst.dsn, self.profdir, 'PNFP', pnfp)
                else:
                    if file_search==False:
                        pnfp = {'knots_spm': self.pnfp_splrep[0],
                                'spline_coef': self.pnfp_splrep[1],
                                'spline_order': self.pnfp_splrep[2]}
                        write_intermediate_files(self.__rsr_inst.year,
                                self.__rsr_inst.doy, self.__rsr_inst.band,
                                self.__rsr_inst.dsn, self.profdir, 'PNFP', pnfp)

            else:
                if (pnfp_coef_ori == self.pnfp_splrep[1]) == False:
                    print('\tPower normalization fit changed within GUI!\n'
                            + '\tWriting new file...')
                    pnfp = {'knots_spm': self.pnfp_splrep[0],
                            'spline_coef': self.pnfp_splrep[1],
                            'spline_order': self.pnfp_splrep[2]}
                    write_intermediate_files(self.__rsr_inst.year,
                            self.__rsr_inst.doy, self.__rsr_inst.band,
                            self.__rsr_inst.dsn, self.profdir, 'PNFP', pnfp)



        self.__set_history()

        return spm_fit, spline_fit#, spline_rep

    def __get_default_freespace_ind(self, rho_km_vals_down, freespace_km,
            is_blocked_atm):
        """
        Purpose:
        Define default freespace regions. Overridden if freespace_spm is not
        None, or if there's an error in the input

        Args:
            rho_km_vals_down (np.ndarray):
                Ring radius values downsampled from raw resolution
            freespace_km (list):
                Default ring regions to include in a fit
            is_blocked_atm (np.ndarray):
                Indices of rho_km_vals_down to ignore
                due to atmosphere or ionosphere occultation
        """

        # Define freespace regions. Overridden if freespace_spm is not None
        ind_free = []
        for i in range(len(freespace_km)):
            ind_free.append(np.argwhere((rho_km_vals_down > freespace_km[i][0])
                & (rho_km_vals_down < freespace_km[i][1])
                & (np.invert(is_blocked_atm))))
        ind_free = np.reshape(np.concatenate(ind_free), -1)

        return ind_free

    def __get_ind_from_knots(self, type_of_input, knots_where, knots,
            spm_vals_free, rho_km_vals_free):
        """
        Purpose:
        Get data points closest to selected knots

        Args:
            type_if_input (str):
                Whether knots_where and knots is input in SPM
                or radius. Is always either "SPM" or "RHO_KM"
            knots_where (np.ndarray):
                Indices of "knots" that are within the specified SPM regions
            knots (list):
                Specified knots to use in the spline fit
            spm_vals_free (np.ndarray):
                SPM values in the freespace regions
            rho_km_vals_free (np.ndarray):
                Radius values in the freespace regions
        """

        # Select data points closest to selected knots
        n_knots = len(knots_where)
        knots_spm_data = np.zeros(n_knots)
        knots_rho_data = np.zeros(n_knots)
        for i in range(n_knots):
            _knot = knots[int(knots_where[i])]
            if type_of_input == 'SPM':
                _ind = np.argmin(abs(spm_vals_free - _knot))
            elif type_of_input == 'RHO_KM':
                _ind = np.argmin(abs(rho_km_vals_free - _knot))
            knots_spm_data[i] = spm_vals_free[_ind]
            knots_rho_data[i] = rho_km_vals_free[_ind]

        return knots_spm_data, knots_rho_data

    def __set_history(self):
        """
        Purpose:
        Record info about the run's history
        """

        input_var_dict = {'geo_inst': self.__geo_inst.history,
            'rsr_inst': self.__rsr_inst.history}
        input_kw_dict = {'spm_fit': self._spm_fit,
            'spline_order': self._spline_order,
            'dt_down': self._dt_down,
            'freespace_spm': self._freespace_spm, 'knots_spm': self._knots_spm}
        hist_dict = rss.tools.write_history_dict(
            input_var_dict, input_kw_dict, __file__)
        self.history = hist_dict
