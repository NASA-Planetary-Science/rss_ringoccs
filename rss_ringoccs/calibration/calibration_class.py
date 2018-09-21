"""

calibration_class.py

Purpose: Class for calibration parameters. This doesn't really do anything. It
         just puts everything calibration-related into one place

WHERE TO GET NECESSARY INPUT:
    fit_inst: Use an instance of the FreqOffsetFit class, found inside of
        rss_ringoccs/calibration/freq_offset_fit.py
    norm_inst: Use an instance of the Normalization class, found inside of
        rss_ringoccs/calibration/power_normalization.py
    geo_inst: Use an instance of the Geometry class, found inside of
        rss_ringoccs/occgeo/calc_occ_geometry.py

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
"""


import numpy as np
from scipy.interpolate import interp1d
import sys
import pdb
import pickle

sys.path.append('../..')
import rss_ringoccs as rss
from rss_ringoccs.tools.search_for_file import search_for_file
from rss_ringoccs.tools.write_intermediate_files import write_intermediate_files
from rss_ringoccs.tools.write_output_files import write_output_files

sys.path.remove('../..')


from .freq_offset_fit import FreqOffsetFit
class Calibration(object):
    """
    Purpose:
    Clump together everything calibration-related into a calibration instance.
    Use this if you haven't made a calibration file yet. If you have made a
    calibration file, use the "MakeCalInst" class instead.

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
        >>> cal_inst = rss.calibration.Calibration(fit_inst, norm_inst,
                geo_inst, dt_cal=dt_cal, verbose=verbose)

    Attributes:
        t_oet_spm_vals (np.ndarray):
            SPM values from calibration file
        f_sky_hz_vals (np.ndarray):
            Predicted sky frequency values from
            calibration file, plus the fit to frequency offset
        f_sky_resid_fit_vals (np.ndarray):
            Fit to residual frequency from calibration file
        p_free_vals (np.ndarray):
            Freespace power spline fit values from calibration file
        f_offset_fit_vals (np.ndarray):
            Fit to frequency offset from calibration files
        history (dict):
            Dictionary with information of the run
    """

    def __init__(self, rsr_inst, geo_inst, dt_cal=1.0, file_search=True,
            verbose=False, USE_GUI=False, write_file=True):
        """
        Args:
            fit_inst:
                Instance of the FreqOffsetFit class. Used for the
                frequency offset fit and the residual frequency fit
            norm_inst:
                Instance of the Normalization class. Used for the power
                normalizing spline fit
            geo_inst:
                Instance of the Geometry class. Used for various
                geometry parameters
            dt_cal (float):
                Desired time Spacing between points
            verbose (bool):
                Print intermediate steps and results

        Dependencies:
            [1] FreqOffsetFit
            [2] Normalization
            [3] Geometry
            [4] scipy
            [5] numpy
        """

     #   if not isinstance(fit_inst, rss.calibration.FreqOffsetFit):
     #       sys.exit('ERROR (Calibration): fit_inst input needs to be an '
     #           + 'instance of the FreqOffsetFit class')

     #   if not isinstance(norm_inst, rss.calibration.Normalization):
     #       sys.exit('ERROR (Calibration): norm_inst input needs to be an '
     #           + 'instance of the Normalization class')

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
                file_search=file_search, USE_GUI=USE_GUI, verbose=verbose)


        # Get corrected I's and Q's
        spm_vals, IQ_c = fit_inst.get_IQ_c()



        spm_geo = np.asarray(geo_inst.t_oet_spm_vals)
        rho_km_geo = np.asarray(geo_inst.rho_km_vals)
        #phi_rad_vals = np.radians(np.asarray(geo_inst.phi_ora_deg_vals))
        #B_rad_vals = np.radians(np.asarray(geo_inst.B_deg_vals))
        #D_km_vals = np.asarray(geo_inst.D_km_vals)
        #rho_dot_kms_vals = np.asarray(geo_inst.rho_dot_kms_vals)
        #F_km_vals = np.asarray(geo_inst.F_km_vals)

        if verbose:
            print('\tRetrieving predicted sky frequency, residual '
                + 'frequency fit, and frequency offset fit...')
        f_spm, f_sky_pred = fit_inst.get_f_sky_pred()
        f_spm, f_sky_resid_fit = fit_inst.get_f_sky_resid_fit()
        f_spm, f_offset_fit = fit_inst.get_f_offset_fit()

        # SPM for calibration parameters
        if verbose:
            print('\tCreating set of SPM at time spacing ' + str(dt_cal)
                    + '...')
        n_pts_cal = round((spm_geo[-1] - spm_geo[0]) / dt_cal) + 1
        spm_cal = spm_geo[0] + dt_cal * np.arange(n_pts_cal)

        if verbose:
            print('\tInterpolating values to calibration SPM...')
                #predicted sky frequency, residual sky '
                #+ 'frequency fit, and frequency offset fit to calibration '
                #+ 'SPM. Evaluating spline fit at calibration SPM')

        # Evaluate f_sky_pred at spm_cal
        f_sky_pred_func = interp1d(f_spm, f_sky_pred, fill_value='extrapolate')
        f_sky_pred_cal = f_sky_pred_func(spm_cal)

        # Evaluate f_sky_resid_fit at spm_cal
        f_sky_resid_fit_func = interp1d(f_spm, f_sky_resid_fit,
            fill_value='extrapolate')
        f_sky_resid_fit_cal = f_sky_resid_fit_func(spm_cal)

        # Evaluate f_offset_fit at spm_cal
        f_offset_fit_func = interp1d(f_spm, f_offset_fit,
            fill_value='extrapolate')
        f_offset_fit_cal = f_offset_fit_func(spm_cal)

        # Normalize observed power by the freespace signal
        norm_inst = rss.calibration.Normalization(spm_vals, IQ_c,
            geo_inst, rsr_inst, verbose=verbose, file_search=file_search)

        profdir = geo_inst.get_profile_dir()
        # If set, search for power normalization fit pickle (PNFP) file
        #if file_search:
        #    pnfp_file = search_for_file(rsr_inst.year,
        #            rsr_inst.doy, rsr_inst.band, 
        #            (rsr_inst.dsn).split('-')[-1], profdir, 'PNFP')

        #    print('\tExtracting power normalization fit from:\n\t\t'
        #            + '/'.join(pnfp_file.split('/')[0:5]) + '/\n\t\t\t'
        #            + pnfp_file.split('/')[-1])
        #    file_object = open(pnfp_file, 'rb')
        #    fit_param_dict = pickle.load(file_object)
        #    k_power_norm = fit_param_dict['k']
        #    freespace_spm = fit_param_dict['freespace_spm']
        #    knots_spm = fit_param_dict['knots_spm']
        #    
        #    spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
        #        freespace_spm=freespace_spm, knots_spm=knots_spm,
        #        spline_order=k_power_norm, USE_GUI=False, verbose=verbose)
            

            

        # Evaluate spline fit at spm_cal. Assumes you already made a
        #     satisfactory spline fit
        dummy_spm, _p_free = norm_inst.get_spline_fit(USE_GUI=USE_GUI,
                file_search=file_search)
        p_free_cal_func = interp1d(dummy_spm, _p_free)
        p_free_cal = p_free_cal_func(spm_cal)

        self.t_oet_spm_vals = spm_cal
        self.f_sky_hz_vals = f_sky_pred_cal + f_offset_fit_cal
        self.f_sky_resid_fit_vals = f_sky_resid_fit_cal
        self.p_free_vals = p_free_cal
        self.__set_history(fit_inst, norm_inst, geo_inst, dt_cal)

        #if not file_search:
        #    pnfp = {'k': norm_inst._spline_order,
        #            'freespace_spm': norm_inst._freespace_spm,
        #            'knots_spm': norm_inst._knots_spm}
        #    write_intermediate_files(rsr_inst.year, rsr_inst.doy,
        #            rsr_inst.band, rsr_inst.dsn, profdir, 'PNFP', pnfp)

        if write_file:
            write_output_files(self)

    def __set_history(self, fit_inst, norm_inst, geo_inst, dt_cal):
        """
        Purpose:
        Set history attribute of the class

        Args:
            fit_inst:
                Instance of the FreqOffsetFit class. Used for the
                frequency offset fit and the residual frequency fit
            norm_inst:
                Instance of the Normalization class. Used for the power
                normalizing spline fit
            geo_inst:
                Instance of the Geometry class. Used for various
                geometry parameters
            dt_cal (float):
                Desired time Spacing between points
        """

        input_var_dict = {'fit_inst': fit_inst.history,
            'norm_inst': norm_inst.history, 'geo_inst': geo_inst.history}
        input_kw_dict = {'dt_cal': dt_cal}
        hist_dict = rss.tools.write_history_dict(
            input_var_dict, input_kw_dict, __file__)

        self.history = hist_dict
