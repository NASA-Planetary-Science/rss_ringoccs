"""
:Purpose:

    Class framework for performing the necessary calibration steps for the RSR
    data. This includes phase correction based on frequency offset of the
    spacecraft and normalization of received power with respect to the intrinsic
    spacecraft power.

:Notes:

    Can be computationally cumbersome, especially for chord occultations. May
    require up to 30 mins for 16 kHz RSR data files.

:Dependencies:

    #. numpy
    #. pdb
    #. scipy
    #. sys

"""


import numpy as np
from scipy.interpolate import splrep,splev
import sys
import pdb

sys.path.append('../..')
import rss_ringoccs as rss
from rss_ringoccs.tools.search_for_file import search_for_file
from rss_ringoccs.tools.write_output_files import write_output_files
from rss_ringoccs.tools.write_history_dict import write_history_dict
from .calc_tau_thresh import calc_tau_thresh
sys.path.remove('../..')


from .freq_offset_fit import FreqOffsetFit
class Calibration(object):
    """
    Purpose:
        Define a class which, when instantiated, calls the submodules for performing
        each step of the calibration process by instantiating the classes
        ``FreqOffsetFit`` in the ``freq_offset_fit.py`` script and ``Normalization``
        in the ``power_normalization.py`` script.
        .. If the calibration files for the desired input parameters already exsit,
        .. use the ``MakeCalInst`` class in the ``cal_inst_from_file.py``.

    Arguments
        :rsr_inst (*object*): Instance of the RSRReader class
        :geo_inst (*object*): Instance of the Geometry class

    Keyword Arguments
        :fof_order (*float*): Whole number specifying the polynomial order to use
                                when fitting the frequency offset residual.
                                Default is 9.
        :pnf_order (*float*): whole number specifying the polynomial order to use
                                when fitting the freespace power. Default is 3.
        :dt_cal (*float*): Desired final spacing in SPM between data points.
                                Default is 1 sec.
        :verbose (*bool*): If True, print intermediate steps and results. Default
                                is False.
        :write_file (*bool*): If True, write output *CAL.TAB and *CAL.LBL files.
                                Default is True.
        :interact (*bool*): If True, enables the interactive mode in the terminal
                                for fitting the freespace power. Default is False.

    Attributes
        :rev_info (*list*): list of strings identifying the specific occultation:
                                year, day of year, direction and type of occultation,
                                spacecraft revolution number, and observation band
        :t_oet_spm_vals (*np.ndarray*): SPM values for occultation event time :math:`t`
        :f_sky_hz_vals (*np.ndarray*): sum of the predicted sky frequency values
                                        and the fit to frequency offset
                                        :math:`\\hat{f}(t)_{offset}=(f_{dr}(t)-f_{dp}(t))+\\hat{f}(t)_{resid}`
        :f_sky_resid_fit_vals (*np.ndarray*): fit to residual sky frequency :math:`\\hat{f}(t)_{resid}`
        :p_free_vals (*np.ndarray*): fit to freespace power :math:`\\hat{P}_0(t)`
        :IQ_c (*np.ndarray*): phase-corrected spacecraft signal :math:`I_{c}+iQ_{c}`
        :history (*dict*): information about the parameters, results, and
                            computation of the calibration procedures

    Returns
        Object instance with attributes associated with the process of calibrating the
        measured signal :math:`I_{m}+iQ_{m}` as well as the method for
        phase-correcting :math:`I_{m}+iQ_{m}`.
    """

    def __init__(self, rsr_inst, geo_inst, fof_order=9, pnf_order=3, dt_cal=1.0,
                 verbose=False, write_file=True, interact=False, pnf_fittype='poly'):

        if not isinstance(geo_inst, rss.occgeo.Geometry):
            sys.exit('ERROR (Calibration): geo_inst input needs to be an '
                + 'instance of the Geometry class')

        dt_cal = abs(dt_cal)
        if (not isinstance(dt_cal, float)) and (not isinstance(dt_cal, int)):
            print('WARNING (Calibration): dt_cal input must be either a '
                + 'float or an integer. Setting to default of 1.0')
            dt_cal = 1.0

        if not isinstance(verbose, bool):
            print('WARNING (Calibration): verbose input must be one of '
                + 'Python\'s built-in booleans (True or False). Setting to '
                + 'False')
            verbose = False

        if verbose:
            print('\nCalibrating frequency and power...')

        # Extract rev_info from geo_inst -- this will inherit its prof_dir
        self.rev_info = geo_inst.rev_info


        # Calculate frequency offset fit
        ## Use default residual frequency fit
        fit_inst = rss.calibration.FreqOffsetFit(rsr_inst, geo_inst,
                verbose=verbose, poly_order=fof_order)

        # Get corrected I's and Q's
        self.IQ_c = self.correct_IQ(rsr_inst.spm_vals,rsr_inst.IQ_m,
                            fit_inst.f_spm,fit_inst.f_offset_fit)

        # Normalize observed power by the freespace signal
        norm_inst = rss.calibration.Normalization(rsr_inst.spm_vals, self.IQ_c,
            geo_inst, rsr_inst, order=pnf_order,verbose=verbose,interact=interact)

        '''
            Resample calibration results to spm_cal
        '''
        spm_cal = np.arange(geo_inst.t_oet_spm_vals[0],geo_inst.t_oet_spm_vals[-1],dt_cal)
        # Evaluate f_sky_pred at spm_cal
        f_sky_pred_splcoef = splrep(fit_inst.f_spm, fit_inst.f_sky_pred)
        f_sky_pred_cal = splev(spm_cal,f_sky_pred_splcoef)
        # Evaluate f_sky_resid_fit at spm_cal
        f_sky_resid_fit_splcoef = splrep(fit_inst.f_spm,fit_inst.f_sky_resid_fit)
        f_sky_resid_fit_cal = splev(spm_cal,f_sky_resid_fit_splcoef)
        # Evaluate f_offset_fit at spm_cal
        f_offset_fit_splcoef = splrep(fit_inst.f_spm,fit_inst.f_offset_fit)
        f_offset_fit_cal = splev(spm_cal,f_offset_fit_splcoef)
        # Evaluate spline fit at spm_cal
        p_free_cal_splcoef = splrep(norm_inst.spm_down, norm_inst.pnorm_fit)
        #p_free_cal_splcoef = splrep(rsr_inst.spm_vals, norm_inst.pnorm_fit)
        p_free_cal = splev(spm_cal,p_free_cal_splcoef)
        p_free_raw = splev(rsr_inst.spm_vals,p_free_cal_splcoef)

        # attributes for writing to file
        self.t_oet_spm_vals = spm_cal
        self.f_sky_hz_vals = f_sky_pred_cal + f_offset_fit_cal
        self.f_sky_resid_fit_vals = f_sky_resid_fit_cal
        self.p_free_vals = p_free_cal
        gaps_used = norm_inst.gaps

        tau_thresh_inst = calc_tau_thresh(rsr_inst,geo_inst,norm_inst.gaps,p_free_raw,res_km=1)
        self.snr = tau_thresh_inst.snr
        self.tau_thresh = tau_thresh_inst.tau_thresh
        self.spm_thresh = tau_thresh_inst.spm_vals
        self.rho_thresh = tau_thresh_inst.rho_vals

        # Write to file
        input_vars = {
                "rsr_inst": rsr_inst.history,
                "geo_inst": geo_inst.history}
        input_kwds = {
                "fof_order": fof_order,
                "pnf_order": pnf_order,
                "dt_cal": dt_cal,
                "pnf_fittype": pnf_fittype,
                "freespace_spm": gaps_used,
                "interact": interact}

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        self.naif_toolkit_version = geo_inst.naif_toolkit_version


        if write_file:
            write_output_files(self)


    def correct_IQ(self,spm_vals,IQ_m,f_spm,f_offset_fit):
        """
        Purpose:

        Apply frequency offset fit to raw measured signal using the signal
        frequencies calculated by ``FreqOffsetFit``. First resamples the frequency
        offset fit to a 0.1 sec separation. Then, computes detrending function
        by integrating frequency offset fit to get phase detrending function :math:`\\psi`.
        Finally, applies phase detrending correction to signal to raw signal such that

        .. math::

            I_{c}+iQ_{c} = [I_{m}+iQ_{m}] \\exp(-i\\psi)

        Arguments:
            :spm_vals (*np.ndarray*): raw SPM values
            :IQ_m (*np.ndarray*): raw complex signal measured by DSN
            :f_spm (*np.ndarray*): SPM sampled for frequency offset calculation
                                in the ``calc_freq_offset`` class in the
                                ``calc_freq_offset.py`` script.
            :f_offset_fit (*np.ndarray*): frequency of the spacecraft signal
                                corresponding to ``f_spm``

        Returns:
            :IQ_c (*np.ndarray*):
                Frequency-corrected complex signal :math:`I_{c}+iQ_{c}`
                corresponding to ``spm_vals``
        """

        # Interpolate frequeny offset fit to 0.1 second spacing, since
        # this makes the integration later more accurate
        dt = 0.1
        npts = round((f_spm[-1] - f_spm[0]) / dt)
        f_spm_interp = f_spm[0] + dt * np.arange(npts)
        f_offset_fit_splcoef = splrep(f_spm, f_offset_fit)
        f_offset_fit_interp = splev(f_spm_interp,f_offset_fit_splcoef)
        #f_offset_fit_function = interp1d(f_spm, f_offset_fit)
        #f_offset_fit_interp = f_offset_fit_function(f_spm_interp)

        # Integration of frequency offset fit to get phase detrending function.
        # Then interpolated to same SPM as I and Q
        f_detrend_interp = np.cumsum(f_offset_fit_interp) * dt
        f_detrend_interp_rad = f_detrend_interp * (2.0 * np.pi)
        f_detrend_rad_splcoef = splrep(f_spm_interp, f_detrend_interp_rad)
        f_detrend_rad = splev(spm_vals, f_detrend_rad_splcoef)
        #f_detrend_rad_function = interp1d(f_spm_interp, f_detrend_interp_rad,
        #        fill_value='extrapolate')
        #f_detrend_rad = f_detrend_rad_function(spm_vals)

        # Apply detrending function
        IQ_c = IQ_m * np.exp(-1j * f_detrend_rad)

        return IQ_c
"""
Revisions:
        calibration_class.py
    2018 Jun 11 - gsteranka - Original version
    2018 Jun 27 - gsteranka - Adjust so sky frequency as frequency offset fit
                              added on top of predicted sky frequency from
                              RSR file
    2018 Sep 16 - jfong - remove fit_inst, norm_inst inputs (create them here)
    2018 Sep 17 - jfong - add file_search kwd
    2018 Sep 18 - jfong - update verbose print statements
    2018 Sep 20 - jfong - add write_file kwarg
                        - set rev_info from geo_inst as attribute
                            - this will inherit the same prof_dir
    2018 Sep 25 - jfong - allow file_search to be a boolean or a string
                            - if string, parse for FRFP or PNFP file paths,
                            - if only one file path given, search for the most
                              recent file for the other file
                        - file_search=False default in __init__ kwarg
    2018 Oct 11 - sflury - update to use only splrep,splev for interpolating
                           and updated to match inner workings of fit_inst
                           and norm_inst -- lots of simplying, Thoreau would be
                           happy
    2018 Nov 15 - sflury - removed "order" kwarg and added "fof_order", "pnf_order",
                           and "interact" kwargs for customizing fits when creating
                           instance of Calibration class with optional interactive
                           mode
    2018 Nov 26 - sflury - reformatted doc strings to match new, standardized format
    2018 Nov 29 - sflury - added call to threshold optical depth class
"""
