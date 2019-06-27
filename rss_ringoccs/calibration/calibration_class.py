"""

:Purpose:
    Class framework for performing the necessary calibration steps for
    the RSR data. This includes phase correction based on frequency
    offset of the spacecraft and normalization of received power with
    respect to the intrinsic spacecraft power.

:Notes:
    Can be computationally cumbersome, especially for chord
    occultations. May require up to 30 mins for 16 kHz RSR data files.

Dependencies:
    #. numpy
    #. scipy
    #. sys

"""

import numpy as np
from scipy.interpolate import splrep,splev
from scipy.integrate import simps
import sys

from ..tools.history import write_history_dict
from ..tools.write_output_files import write_output_files
from .freq_offset_fit import FreqOffsetFit
from .power_normalization import Normalization
from ..occgeo import Geometry

class Calibration(object):
    """
    :Purpose:
        Define a class which, when instantiated, calls the submodules
        for performing each step of the calibration process by
        instantiating the classes ``FreqOffsetFit`` in the
        ``freq_offset_fit.py`` script and ``Normalization`` in the
        ``power_normalization.py`` script.

    Arguments
        :rsr_inst (*object*): Instance of the RSRReader class
        :geo_inst (*object*): Instance of the Geometry class

    Keyword Arguments
        :pnf_order (*float*): whole number specifying the polynomial
                        order to use when fitting the freespace power.
                        Default is 3.
        :dt_cal (*float*): Desired final spacing in SPM between data
                        points. Default is 1 sec.
        :verbose (*bool*): If True, print intermediate steps and
                        results. Default is False.
        :write_file (*bool*): If True, write output CAL .TAB and
                        CAL .LBL files. Default is True.
        :interact (*bool*): If True, enables the interactive mode in
                        the terminal for fitting the freespace power.
                        Default is False.

    Attributes
        :rev_info (*dict*): *dict* of information identifying the
                        specific occultation: rsrfile, year, day of
                        year, direction and type of occultation,
                        spacecraft revolution number, and observation
                        band
        :t_oet_spm_vals (*np.ndarray*): SPM values for observed event
                        time :math:`t`
        :f_sky_hz_vals (*np.ndarray*):
                        sum of the reconstructed sky frequency values
                        and the fit to frequency offset
                        :math:`\\hat{f}(t)_{offset}=(f_{dr}(t)
                        -f_{dp}(t))+\\hat{f}(t)_{offset}`
                        following Equation 19 in [CRSUG2018]_.
        :f_offset_fit_vals (*np.ndarray*): fit to
                        frequency offset :math:`\\hat{f}(t)_{offset}`
        :p_free_vals (*np.ndarray*): fit to freespace power
                        :math:`\\hat{P}_0(t)`
        :IQ_c (*np.ndarray*): phase-corrected spacecraft signal
                        :math:`I_{c}+iQ_{c}`
        :history (*dict*): information about the parameters, results,
                        and computation of the calibration procedures
        :FORFIT_chi_squared (*float*): sum of the squared residual
                        frequency offset fit such that
                        :math:`\chi^2 = \\frac{1}{N-m}
                        \sum((\hat{f}(t)_{offset}-f(t)_{offset})
                        /\hat{f}(t)_{offset})^2`
        :FSPFIT_chi_squared (*float*):
                        :math:`\chi^2 = \\frac{1}{N-m}\sum
                        ((\hat{P}_0(t)-P_0(t))/\hat{P}_0(t))^2`
    """

    def __init__(self, rsr_inst, geo_inst, pnf_order=3, dt_cal=1.0,
                 verbose=False, write_file=True, interact=False):

        if not isinstance(geo_inst, Geometry):
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

        '''
        # Find the max and min SPM values with True boolean indices
        if len(self.raw_spm_vals[inds]) > 2 :
            # if user-defined lower spm is more constraining than geometry
            if np.min(self.raw_spm_vals[inds]) < np.min(spm_limits) :
                spm_min = np.min(spm_limits)
            # otherwise, use limits implied by geometry
            else:
                spm_min = np.min(self.raw_spm_vals[inds])
            # if user-defined upper spm is more constraining than geometry
            if np.max(self.raw_spm_vals[inds]) > np.max(spm_limits) :
                spm_max = np.max(spm_limits)
            # otherwise, use limits implied by geometry
            else:
                spm_max = np.max(self.raw_spm_vals[inds])
        # if geometry is too restricting, use max and min SPM
        else:
            print('Error in estimating SPM range for frequency offset!')
            spm_min = self.raw_spm_vals[0]
            spm_max = self.raw_spm_vals[-1]
        '''

        if verbose:
            print('\nCalibrating frequency and power...')

        # Extract rev_info from geo_inst -- this will inherit its prof_dir
        self.rev_info = geo_inst.rev_info


        # Calculate frequency offset fit
        # Use default offset frequency fit
        fit_inst = FreqOffsetFit(rsr_inst, geo_inst, verbose=verbose,
                write_file=write_file)

        # Get corrected I's and Q's
        self.IQ_c = self.correct_IQ(rsr_inst.spm_vals,rsr_inst.IQ_m,
                            fit_inst.f_spm,fit_inst.f_offset_fit)

        # Normalize observed power by the freespace signal
        norm_inst = Normalization(rsr_inst.spm_vals, self.IQ_c, geo_inst,
                order=pnf_order,verbose=verbose,interact=interact,
                write_file=write_file)

        spm_cal = np.arange(geo_inst.t_oet_spm_vals[0],
                geo_inst.t_oet_spm_vals[-1],dt_cal)

        # Evaluate f_sky_recon at spm_cal
        f_sky_recon_splcoef = splrep(fit_inst.f_spm, fit_inst.f_sky_recon)
        f_sky_recon_cal = splev(spm_cal,f_sky_recon_splcoef)

        # Evaluate f_offset_fit at spm_cal
        f_offset_fit_splcoef = splrep(fit_inst.f_spm,fit_inst.f_offset_fit)
        f_offset_fit_cal = splev(spm_cal,f_offset_fit_splcoef)

        # Evaluate spline fit at spm_cal
        p_free_cal_splcoef = splrep(norm_inst.spm_down, norm_inst.pnorm_fit)
        p_free_cal = splev(spm_cal,p_free_cal_splcoef)
        p_free_raw = splev(rsr_inst.spm_vals,p_free_cal_splcoef)

        # attributes for writing to file
        self.t_oet_spm_vals = spm_cal
        self.f_sky_hz_vals = f_sky_recon_cal + f_offset_fit_cal
        self.f_offset_fit_vals = f_offset_fit_cal
        self.f_offset = fit_inst.f_offset
        self.f_spm = fit_inst.f_spm
        self.p_free_vals = p_free_cal
        gaps_used = norm_inst.gaps
        self.gaps = norm_inst.gaps
        self.FORFIT_chi_squared = fit_inst.chi_squared
        self.FSPFIT_chi_squared = norm_inst.chi_squared

        # Write to file
        input_vars = {
                "rsr_inst": rsr_inst.history,
                "geo_inst": geo_inst.history}
        input_kwds = {
                "pnf_order": pnf_order,
                "dt_cal": dt_cal,
                "interact": interact}

        additional_info = {
                "FORFIT_chi_squared": self.FORFIT_chi_squared,
                "FSPFIT_chi_squared": self.FSPFIT_chi_squared,
                "pnf_fittype": 'poly',
                "fof_order": fit_inst.poly_order,
                "freespace_spm": gaps_used}


        self.history = write_history_dict(input_vars, input_kwds, __file__,
                add_info=additional_info)
        self.naif_toolkit_version = geo_inst.naif_toolkit_version


        if write_file:
            self.outfiles = write_output_files(self)


    def correct_IQ(self,spm_vals,IQ_m,f_spm,f_offset_fit):
        """
        Purpose:

            Apply frequency offset fit to raw measured signal using
            the signal frequencies calculated by ``FreqOffsetFit``.
            First resamples the frequency offset fit to a 0.1 sec
            separation. Then, computes detrending function by
            integrating frequency offset fit to get phase detrending
            function :math:`\\psi` using Equation 18 from
            [CRSUG2018]_ where

            .. math::
                \\psi = \int^{t}\hat{f}(\\tau)_{offset}
                \mathrm{d}\\tau+\\psi(t_0)

            Finally, applies phase detrending correction to signal to
            raw signal such that

            .. math::
                I_{c}+iQ_{c} = [I_{m}+iQ_{m}] \\exp(-i\\psi)

            as discussed in [CRSUG2018]_ (see their Equation 17).

        Arguments:
            :spm_vals (*np.ndarray*): raw SPM values
            :IQ_m (*np.ndarray*): raw complex signal measured by DSN
            :f_spm (*np.ndarray*): SPM sampled for frequency offset
                        calculation in the ``calc_freq_offset`` class
                        in the ``calc_freq_offset.py`` script.
            :f_offset_fit (*np.ndarray*): frequency of the spacecraft
                        signal corresponding to ``f_spm``

        Returns:
            :IQ_c (*np.ndarray*): Frequency-corrected complex signal
                        :math:`I_{c}+iQ_{c}` corresponding to
                        ``spm_vals``
        """

        # Interpolate frequeny offset fit to 0.1 second spacing, since
        # this makes the integration later more accurate
        dt = 0.01
        npts = round((f_spm[-1] - f_spm[0]) / dt)
        f_spm_ip = f_spm[0] + dt * np.arange(npts)
        f_off_splcoef = splrep(f_spm, f_offset_fit)
        f_off_ip = splev(f_spm_ip,f_off_splcoef)

        # Integration of frequency offset fit to get phase detrending function.
        f_detrend_ip = np.cumsum(f_off_ip)*dt
        #f_detrend_ip = np.zeros(len(f_off_ip))
        #for i in range(len(f_off_ip)):
        #    f_detrend_ip[i] = simps(f_off_ip[:i+1],x=f_spm_ip[:i+1])
        f_detrend_ip_rad = f_detrend_ip * (2.0 * np.pi)
        # Then interpolate to same SPM as I and Q
        f_detrend_rad_splcoef = splrep(f_spm_ip, f_detrend_ip_rad)
        f_detrend_rad = splev(spm_vals, f_detrend_rad_splcoef)

        # Apply detrending function
        IQ_c = IQ_m * np.exp(-1j * f_detrend_rad)

        return IQ_c
"""
Revisions:

"""
