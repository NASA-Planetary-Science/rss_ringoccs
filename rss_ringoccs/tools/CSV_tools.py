"""
    Purpose:
        Provide tools for reading in .TAB and
        .CSV files and converting the data into
        a usable instance of the DLP class.
    Dependencies:
        #. pandas
        #. numpy
        #. scipy
        #. rss_ringoccs
"""

import numpy as np
import pandas as pd
from scipy import interpolate
from .history import write_history_dict, date_to_rev, rev_to_occ_info
from . import error_check

def get_geo(geo, verbose=True, use_deprecate=False):
    """
        Purpose:
            To extract a pandas DataFrame from a given
            GEO.TAB or GEO.CSV file.
        Arguments:
            :geo (*str*):
                A string containing the location of
                the requested geo file.
                Ex: geo = "/path/to/geo.CSV"
                File must contain the following columns,
                in the following order:
                |   t_oet_spm_vals
                |   t_ret_spm_vals
                |   t_set_spm_vals
                |   rho_km_vals
                |   phi_rl_deg_vals
                |   phi_ora_deg_vals
                |   B_deg_vals
                |   D_km_vals
                |   rho_dot_kms_vals
                |   phi_rl_dot_kms_vals
                |   F_km_vals
                |   R_imp_km_vals
                |   rx_km_vals
                |   ry_km_vals
                |   rz_km_vals
                |   vx_kms_vals
                |   vy_kms_vals
                |   vz_kms_vals
                |   obs_spacecract_lat_deg_vals
        Keywords:
            :verbose (*bool*):
                A Boolean for printing out auxiliary
                information to the command line.
    """
    fname = "tools.CSV_tools.get_geo"
    error_check.check_type(geo, str, "geo", fname)
    error_check.check_type(verbose, bool, "verbose", fname)

    if verbose:
        print("\tExtracting Geo Data...")

    try:
        if use_deprecate:
            dfg = pd.read_csv(geo, delimiter=',',
                              names=["t_oet_spm_vals",
                                     "t_ret_spm_vals",
                                     "t_set_spm_vals",
                                     "rho_km_vals",
                                     "phi_rl_deg_vals",
                                     "phi_ora_deg_vals",
                                     "B_deg_vals",
                                     "D_km_vals",
                                     "rho_dot_kms_vals",
                                     "phi_rl_dot_kms_vals",
                                     "F_km_vals",
                                     "R_imp_km_vals",
                                     "rx_km_vals",
                                     "ry_km_vals",
                                     "rz_km_vals",
                                     "vx_kms_vals",
                                     "vy_kms_vals",
                                     "vz_kms_vals"])
        else:
            dfg = pd.read_csv(geo, delimiter=',',
                              names=["t_oet_spm_vals",
                                     "t_ret_spm_vals",
                                     "t_set_spm_vals",
                                     "rho_km_vals",
                                     "phi_rl_deg_vals",
                                     "phi_ora_deg_vals",
                                     "B_deg_vals",
                                     "D_km_vals",
                                     "rho_dot_kms_vals",
                                     "phi_rl_dot_kms_vals",
                                     "F_km_vals",
                                     "R_imp_km_vals",
                                     "rx_km_vals",
                                     "ry_km_vals",
                                     "rz_km_vals",
                                     "vx_kms_vals",
                                     "vy_kms_vals",
                                     "vz_kms_vals",
                                     "obs_spacecract_lat_deg_vals"])
    except FileNotFoundError:
        raise FileNotFoundError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                Your input file does not exist.
                Your file: %s
            """ % (fname, geo)
        )

    if verbose:
        print("\tGeo Data Complete.")

    return dfg

def get_cal(cal, verbose=True):
    """
        Purpose:
            To extract a pandas DataFrame from a given
            CAL.TAB or CAL.CSV file.
        Arguments:
            :cal (*str*):
                A string containing the location of
                the requested cal file.
                Ex: cal = "/path/to/cal.CSV"
                File must contain the following columns,
                in the following order:
                |   spm_vals
                |   f_sky_pred_vals
                |   f_sky_resid_fit_vals
                |   p_free_vals
        Keywords:
            :verbose (*bool*):
                A Boolean for printing out auxiliary
                information to the command line.
    """
    fname = "tools.CSV_tools.get_cal"
    error_check.check_type(cal, str, "cal", fname)
    error_check.check_type(verbose, bool, "verbose", fname)

    if verbose:
        print("\tExtracting Cal Data...")

    try:
        dfc = pd.read_csv(cal, delimiter=',',
                          names=["spm_vals",
                                 "f_sky_pred_vals",
                                 "f_sky_resid_fit_vals",
                                 "p_free_vals"])
    except FileNotFoundError:
        raise FileNotFoundError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                Your input file does not exist.
                Your file: %s
            """ % (fname, cal)
        )

    if verbose:
        print("\tCal Data Complete.")

    return dfc

def get_dlp(dlp, verbose=True, use_deprecate=False):
    """
        Purpose:
            To extract a pandas DataFrame from a given
            DLP.TAB or DLP.CSV file.
        Arguments:
            :dlp (*str*):
                A string containing the location of
                the requested dlp file.
                Ex: dlp = "/path/to/dlp.CSV"
                File must contain the following columns,
                in the following order:
                |   rho_km_vals
                |   rho_corr_pole_km_vals
                |   rho_corr_timing_km_vals
                |   phi_rl_deg_vals
                |   phi_ora_deg_vals
                |   p_norm_vals
                |   raw_tau_vals
                |   phase_deg_vals
                |   raw_tau_threshold_vals
                |   t_oet_spm_vals
                |   t_ret_spm_vals
                |   t_set_spm_vals
                |   B_deg_vals
        Keywords:
            :verbose (*bool*):
                A Boolean for printing out auxiliary
                information to the command line.
    """
    fname = "tools.CSV_tools.get_dlp"
    error_check.check_type(dlp, str, "dlp", fname)
    error_check.check_type(verbose, bool, "verbose", fname)

    if verbose:
        print("\tExtracting DLP Data...")

    try:
        if use_deprecate:
            dfd = pd.read_csv(dlp, delimiter=',',
                              names=["rho_km_vals",
                                     "rho_corr_pole_km_vals",
                                     "rho_corr_timing_km_vals",
                                     "phi_rl_deg_vals",
                                     "phi_ora_deg_vals",
                                     "raw_tau_vals",
                                     "phase_deg_vals",
                                     "raw_tau_threshold_vals",
                                     "t_oet_spm_vals",
                                     "t_ret_spm_vals",
                                     "t_set_spm_vals",
                                     "B_deg_vals"])
        else:
            dfd = pd.read_csv(dlp, delimiter=',',
                              names=["rho_km_vals",
                                     "rho_corr_pole_km_vals",
                                     "rho_corr_timing_km_vals",
                                     "phi_rl_deg_vals",
                                     "phi_ora_deg_vals",
                                     "p_norm_vals",
                                     "raw_tau_vals",
                                     "phase_deg_vals",
                                     "raw_tau_threshold_vals",
                                     "t_oet_spm_vals",
                                     "t_ret_spm_vals",
                                     "t_set_spm_vals",
                                     "B_deg_vals"])
    except FileNotFoundError:
        raise FileNotFoundError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                Your input file does not exist.
                Your file: %s
            """ % (fname, dlp)
        )

    if verbose:
        print("\tDLP Data Complete")
    return dfd

def get_tau(tau, verbose=True, use_deprecate=False):
    """
        Purpose:
            To extract a pandas DataFrame from a given
            TAU.TAB or TAU.CSV file.
        Arguments:
            :tau (*str*):
                A string containing the location of
                the requested tau file.
                Ex: tau = "/path/to/tau.CSV"
                File must contain the following columns,
                in the following order:
                |   rho_km_vals
                |   rho_corr_pole_km_vals
                |   rho_corr_timing_km_vals
                |   phi_rl_deg_vals
                |   phi_ora_deg_vals
                |   p_norm_vals
                |   raw_tau_vals
                |   phase_deg_vals
                |   raw_tau_threshold_vals
                |   t_oet_spm_vals
                |   t_ret_spm_vals
                |   t_set_spm_vals
                |   B_deg_vals
        Keywords:
            :verbose (*bool*):
                A Boolean for printing out auxiliary
                information to the command line.
    """
    fname = "tools.CSV_tools.get_tau"
    error_check.check_type(tau, str, "tau", fname)
    error_check.check_type(verbose, bool, "verbose", fname)

    if verbose:
        print("\tExtracting Tau Data...")

    try:
        if use_deprecate:
            dft = pd.read_csv(tau, delimiter=',',
                              names=["rho_km_vals",
                                     "rho_km_pole_corr_vals",
                                     "rho_km_offsett_vals",
                                     "phi_rl_deg_vals",
                                     "phi_ora_deg_vals",
                                     "raw_tau_vals",
                                     "phase_deg_vals",
                                     "raw_tau_threshold_vals",
                                     "spm_vals",
                                     "t_ret_spm_vals",
                                     "t_set_spm_vals",
                                     "B_deg_vals"])
        else:
            dft = pd.read_csv(tau, delimiter=',',
                              names=["rho_km_vals",
                              "rho_km_pole_corr_vals",
                              "rho_km_offsett_vals",
                              "phi_rl_deg_vals",
                              "phi_ora_deg_vals",
                              "p_norm_vals",
                              "raw_tau_vals",
                              "phase_deg_vals",
                              "raw_tau_threshold_vals",
                              "spm_vals",
                              "t_ret_spm_vals",
                              "t_set_spm_vals",
                              "B_deg_vals"])
    except FileNotFoundError:
        raise FileNotFoundError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                Your input file does not exist.
                Your file: %s
            """ % (fname, tau)
        )

    if verbose:
        print("\tTau Data Complete")

    return dft


class ExtractCSVData(object):
    """
        Purpose:
            Read three csv files (Geo, Cal, and DLP) and return
            an instance containing all necessary attributes to run
            diffraction correction. This instance can be fed
            directly into the DiffractionCorrection class.
        Variables:
            :geo (*str*):
                A string that contains the location of
                the requested Geo file.
                Ex: geo = "/path/to/geo.CSV"
            :cal (*str*):
                    A string that contains the location of
                    the requested Cal file.
                    Ex: cal = "/path/to/cal.CSV"
            dlp (*str*):
                A string that contains the location of
                the requested dlp file.
                Ex: dlp = "/path/to/dlp.CSV"
        Keywords:
            :tau (*str*):
                A string that contains the location
                of the requested Tau file. If not set,
                variables from the tau file will have
                NoneType. Ex: tau = "/path/to/tau.CSV"
            :verbose (*bool*):
                A Boolean for specifying if various
                status updates will be printed to the
                command line.
        Attributes:
            :B_rad_vals:
                The ring opening angle of the ring plane
                with respect to the line of sight from
                Earth to the spacecraft.
            :D_km_vals:
                The distance from the spacecraft to the
                ring-intercept point.
            :f_sky_hz_vals:
                The sky frequency of the incoming signal.
            :p_norm_vals:
                The normalized diffracted power.
            :phase_rad_vals:
                The diffracted phase, in radians.
            :phase_vals:
                The reconstructed phase contained in the
                tau file. If tau is not set, this will be
                a NoneType variable. Units are in radians.
            :phi_rad_vals:
                The observed ring azimuth angle, in radians.
            :phi_rl_rad_vals:
                The observed ring longitude angle, in radians.
            :power_vals:
                The reconstructed power contained in the tau
                file. If tau is not set, this will be a
                NoneType variable. Power is normalized to one
                in free space regions.
            :raw_tau_threshold_vals:
                The threshold optical depth corresponding
                to the diffracted optical depth profile.
            :rho_corr_pole_km_vals:
                Corrections for the ring-intercept point
                computed by taking into account
                Saturn's pole direction.
            :rho_corr_timing_km_vals:
                Timing offset corrections to the
                ring-intercept point.
            :rho_dot_kms_vals:
                The rate of change of the ring-intercept
                point as a function of time. That is,
                drho/dt.
            :rho_km_vals:
                The ring intercept point, in kilometers.
            :t_oet_spm_vals:
                Observed event time, the time the signal
                is recieved on Earth, computed in
                Seconds Past Midnight.
            :t_ret_spm_vals:
                Ring event time, the time the signal crosses
                the rings, computed in Seconds
                Past Midnight.
            :t_set_spm_vals:
                Spacecraft Event Time, the time the signal
                was transmitted from the spacecraft,
                computed in Seconds Past Midnight.
            :tau_rho:
                The ring-intercept point corresponding
                to the values in the tau file. If tau
                is not set, this will be a NoneType
                variable. Units are in kilometers.
            :tau_vals:
                The normalized optical depth contained
                in the tau file. If tau is not set, this
                will be a NoneType variable.
    """
    def __init__(self, geo, cal, dlp, tau=None, verbose=True, use_deprecate=False):
        fname = "tools.CSV_tools.ExtractCSVData"
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(cal, str, "cal", fname)
        error_check.check_type(dlp, str, "dlp", fname)
        error_check.check_type(verbose, bool, "verbose", fname)

        if (not isinstance(tau, type(None))):
            error_check.check_type(tau, str, "tau", fname)
        else:
            pass

        if verbose:
            print("Extracting Data from CSV Files:")

        # Save inputs as attributes.
        self.geo = geo
        self.cal = cal
        self.dlp = dlp
        self.tau = tau

        # Extract GEO, CAL, and DLP data.
        geo_dat = get_geo(self.geo, verbose=verbose, use_deprecate=use_deprecate)
        cal_dat = get_cal(self.cal, verbose=verbose)
        dlp_dat = get_dlp(self.dlp, verbose=verbose, use_deprecate=use_deprecate)

        if verbose:
            print("\tRetrieving Variables...")

        try:
            # Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = np.array(dlp_dat.rho_km_vals)
            errmess = "raw_tau_vals"
            self.raw_tau_vals = np.array(dlp_dat.raw_tau_vals)
            errmess = "phase_rad_vals"
            self.phase_rad_vals = np.array(dlp_dat.phase_deg_vals)
            self.phase_rad_vals = np.deg2rad(self.phase_rad_vals)
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = np.array(geo_dat.phi_ora_deg_vals)
            self.phi_rad_vals = np.deg2rad(self.phi_rad_vals)
            errmess = "B_rad_vals"
            self.B_rad_vals = np.array(geo_dat.B_deg_vals)
            self.B_rad_vals = np.deg2rad(self.B_rad_vals)
            errmess = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(dlp_dat.t_oet_spm_vals)
            errmess = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(dlp_dat.t_ret_spm_vals)
            errmess = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(dlp_dat.t_set_spm_vals)
            errmess = "rho_corr_pole_km_vals"
            self.rho_corr_pole_km_vals = np.array(dlp_dat.rho_corr_pole_km_vals)
            errmess = "rho_corr_pole_km_vals"
            self.rho_corr_timing_km_vals = np.array(dlp_dat.rho_corr_timing_km_vals)
            errmess = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = np.array(geo_dat.phi_rl_deg_vals)
            self.phi_rl_rad_vals = np.deg2rad(self.phi_rl_rad_vals)
            errmess = "raw_tau_threshold_vals"
            self.raw_tau_threshold_vals = np.array(dlp_dat.raw_tau_threshold_vals)
            errmess = "geo_rho"
            geo_rho = np.array(geo_dat.rho_km_vals)
            errmess = "D_km_vals"
            self.D_km_vals = np.array(geo_dat.D_km_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(geo_dat.rho_dot_kms_vals)
            errmess = "f_sky_pred"
            f_sky_pred = np.array(cal_dat.f_sky_pred_vals)
            errmess = "f_sky_resid"
            f_sky_resid = np.array(cal_dat.f_sky_resid_fit_vals)
            errmess = "f_sky_resid"
            self.f_sky_hz_vals = f_sky_pred - f_sky_resid
        except (ValueError, TypeError, NameError, AttributeError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not convert %s into a numpy array.
                    \r\tCheck input GEO, CAL, and DLP files.
                """ % (fname, errmess)
            )

        error_check.check_is_real(self.rho_corr_timing_km_vals,
                                  "rho_corr_timing_km_vals", fname)
        error_check.check_is_real(self.raw_tau_threshold_vals,
                                  "raw_tau_threshold_vals", fname)
        error_check.check_is_real(self.rho_corr_pole_km_vals,
                                  "rho_corr_pole_km_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals,
                                  "rho_dot_kms_vals",fname)
        error_check.check_is_real(self.phi_rl_rad_vals,
                                  "phi_rl_rad_vals", fname)
        error_check.check_is_real(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_is_real(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_is_real(self.t_set_spm_vals, "t_set_spm_vals", fname)
        error_check.check_is_real(self.phase_rad_vals, "phase_rad_vals", fname)
        error_check.check_is_real(self.f_sky_hz_vals, "f_sky_hz_vals", fname)
        error_check.check_is_real(self.raw_tau_vals, "raw_tau_vals", fname)
        error_check.check_is_real(self.phi_rad_vals, "phi_rad_vals", fname)
        error_check.check_is_real(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_is_real(self.B_rad_vals, "B_rad_vals", fname)
        error_check.check_is_real(self.D_km_vals, "D_km_vals", fname)
        error_check.check_is_real(geo_rho, "geo_rho", fname)

        error_check.check_positive(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_positive(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_positive(self.t_set_spm_vals, "t_set_spm_vals", fname)
        error_check.check_positive(self.f_sky_hz_vals, "f_sky_hz_vals", fname)
        error_check.check_positive(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_positive(self.D_km_vals, "D_km_vals", fname)
        error_check.check_positive(geo_rho, "geo_rho", fname)

        error_check.check_two_pi(self.phi_rl_rad_vals, "phi_rl_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phase_rad_vals, "phase_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rad_vals, "phi_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.B_rad_vals, "B_rad_vals",
                                 fname, deg=False)

        if verbose:
            print("\tComputing Variables...")

        if (np.size(self.rho_km_vals) != np.size(self.t_oet_spm_vals)):
            raise ValueError("len(rho_km_vals) != len(t_oet_spm_vals")

        dr = np.zeros(np.size(self.rho_km_vals) - 1)
        dt = np.zeros(np.size(self.t_oet_spm_vals) - 1)

        for i in range(np.size(self.rho_km_vals) - 1):
            dr[i] = self.rho_km_vals[i+1] - self.rho_km_vals[i]
            dt[i] = self.t_oet_spm_vals[i+1] - self.t_oet_spm_vals[i]
        
        drdt = dr/dt

        if (np.min(drdt) < 0.0) and (np.max (drdt) > 0.0):
           raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tdrho/dt has positive and negative values.
                    \r\tCheck your DLP file for errors.
                    \r\tYour file: %s
                """ % (fname, dlp)
            )
        elif (np.size((drdt == 0).nonzero()) != 0):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tdrho/dt has zero valued elements.
                    \r\tCheck your DLP file for errors.
                    \r\tYour file: %s
                """ % (fname, dlp)
            )
        elif (drdt < 0.0).all():
            occ = 'ingress'
        elif (drdt > 0.0).all():
            occ = 'egress'
        else:
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not determine what type of occultation
                    \r\tthis is. Set occ='ingress' or occ='egress'.
                """ % (fname)
            )

        if (occ == 'ingress'):
            crange = (self.rho_dot_kms_vals < 0.0).nonzero()
        elif (occ == 'egress'):
            crange = (self.rho_dot_kms_vals > 0.0).nonzero()
        else:
            crange_e = (self.rho_dot_kms_vals > 0.0).nonzero()
            crange_i = (self.rho_dot_kms_vals < 0.0).nonzero()
            n_e = np.size(crange_e)
            n_i = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals has positive and negative values.
                        \r\tThis is likely a chord occultation.
                        \r\tSet occ='ingress' or occ='egress'
                    """ % (fname)
                )
            elif (n_e == 0) and (n_i == 0):
                raise ValueError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals is either zero or empty.
                        \r\tCheck input geo file for errors.
                        \r\tYour file: %s
                    """ % (fname, geo)
                )
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else:
                raise ValueError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\tCould not determine what type of occultation
                        \r\tthis is. Set occ='ingress' or occ='egress'.
                    """ % (fname)
                )

            del n_e, n_i, crange_e, crange_i

        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                raise TypeError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals is never negative.
                    """ % fname
                )
            elif (occ == 'egress'):
                raise TypeError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals is never negative.
                    """ % fname
                )
            else:
                raise TypeError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\tCould not determine occultation type.
                        \r\tCheck your input GEO file.
                    """ % fname
                )
        else:
            pass


        if verbose:
            print("\tInterpolating Data...")

        geo_rho = geo_rho[crange]
        self.D_km_vals = self.D_km_vals[crange]
        self.B_rad_vals = self.B_rad_vals[crange]
        self.phi_rad_vals = self.phi_rad_vals[crange]
        self.phi_rl_rad_vals = self.phi_rl_rad_vals[crange]
        self.rho_dot_kms_vals = self.rho_dot_kms_vals[crange]

        del crange, geo_dat, cal_dat, dlp_dat

        rfin = int(np.max((np.max(geo_rho)-self.rho_km_vals>=0.0).nonzero()))
        rstart = int(np.min((self.rho_km_vals-np.min(geo_rho)>=0.0).nonzero()))
        n_rho_vals = np.size(self.rho_km_vals)
        n_f_vals = np.size(self.f_sky_hz_vals)
        frange = np.arange(n_f_vals)
        xrange = np.arange(n_rho_vals)*(n_f_vals-1.0)/(n_rho_vals-1.0)
        interp = interpolate.interp1d(geo_rho, self.D_km_vals)
        self.D_km_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.rho_dot_kms_vals)
        self.rho_dot_kms_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rad_vals)
        self.phi_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.B_rad_vals)
        self.B_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rl_rad_vals)
        self.phi_rl_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(frange, self.f_sky_hz_vals)
        self.f_sky_hz_vals = interp(xrange)

        del n_rho_vals, n_f_vals, interp, frange, xrange, geo_rho

        raw_mu = np.sin(np.abs(self.B_rad_vals))
        self.p_norm_vals = np.exp(-self.raw_tau_vals/raw_mu)
        self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals[rstart:rfin+1]
        self.rho_corr_pole_km_vals = self.rho_corr_pole_km_vals[rstart:rfin+1]
        self.t_ret_spm_vals = self.t_ret_spm_vals[rstart:rfin+1]
        self.t_set_spm_vals = self.t_set_spm_vals[rstart:rfin+1]
        self.t_oet_spm_vals = self.t_oet_spm_vals[rstart:rfin+1]
        self.phase_rad_vals = self.phase_rad_vals[rstart:rfin+1]
        self.rho_km_vals = self.rho_km_vals[rstart:rfin+1]
        self.p_norm_vals = self.p_norm_vals[rstart:rfin+1]

        del rstart, rfin, raw_mu

        if (not isinstance(tau, type(None))):
            tau_dat = get_tau(self.tau, verbose=verbose, use_deprecate=use_deprecate)
            tm = np.sin(np.abs(np.deg2rad(tau_dat.B_deg_vals)))
            rmin = np.min(tau_dat.rho_km_vals)
            rmax = np.max(tau_dat.rho_km_vals)
            rfin = int(np.max((rmax-self.rho_km_vals>=0).nonzero()))
            rstart = int(np.min((self.rho_km_vals-rmin>=0).nonzero()))
            self.tau_rho = self.rho_km_vals[rstart:rfin+1]
            interp = interpolate.interp1d(tau_dat.rho_km_vals,
                                          tau_dat.raw_tau_vals, kind='cubic')
            self.tau_vals = interp(self.tau_rho)
            interp = interpolate.interp1d(tau_dat.rho_km_vals,
                                          np.deg2rad(tau_dat.phase_deg_vals),
                                          kind='cubic')
            self.phase_vals = interp(self.tau_rho)
            interp = interpolate.interp1d(tau_dat.rho_km_vals, tm, kind='cubic')
            tm = interp(self.tau_rho)
            self.power_vals = np.exp(-self.tau_vals/tm)
            del tau_dat, tm, rmin, rmax, rfin, rstart, interp
        else:
            self.tau_rho = None
            self.tau_vals = None
            self.phase_vals = None
            self.power_vals = None

        if verbose:
            print("\tData Extraction Complete.")

        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data": self.geo,
            "CAL Data": self.cal,
            "DLP Data": self.dlp
        }

        input_kwds = {"TAU Data": self.tau}

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        var = geo.split("/")[-1]

        var = var.split("_")
        band = '"%s"' % var[3][0]
        year = var[1]
        doy = var[2]
        dsn = "DSS-%s" % (var[3][1:])
        rev_num = date_to_rev(int(year), int(doy))
        occ_dir = rev_to_occ_info(rev_num)
        prof_dir = '"%s"' % var[4]

        self.rev_info = {
            "rsr_file": "UNK",
            "band": band,
            "year": year,
            "doy": doy,
            "dsn": dsn,
            "occ_dir": occ_dir,
            "planetary_occ_flag": occ_dir,
            "rev_num": rev_num,
            "prof_dir": prof_dir
        }

        if verbose:
            print("\tHistory Complete.")

        if verbose:
            print("\tExtract CSV Data Complete.")


class GetUranusData(object):
    def __init__(self,geodata,dlpdata,dx=0.25,occ=None,verbose=False):
        if (not isinstance(geodata,str)):
            raise TypeError("geodata must be a string: '/path/to/geodata'")
        if (not isinstance(dlpdata,str)):
            raise TypeError("dlpdata must be a string: '/path/to/dlpdata'")
        if (not isinstance(dx,float)):
            raise TypeError("dx must be a floating point number")
        if (dx <= 0.0):
            raise ValueEorr("dx must be a positive number")
        if occ:
            if (not isinstance(occ,str)):
                raise TypeError("occ must be a string")
            else:
                occ = occ.replace(" ", "").lower()
                if (occ != 'ingress') and (occ != 'egress'):
                    raise ValueError("occ must be 'egress' of 'ingress'")
                else:
                    pass

        geo_dat = get_geo(geodata,verbose=verbose)
        dlp_dat = pd.read_csv(
            dlpdata, delimiter=',',
            names=[
                "t_oet_spm_vals",
                "p_norm_vals",
                "phase_rad_vals",
                "f_sky_hz_vals"
            ]
        )

        dlp_spm = np.array(dlp_dat.t_oet_spm_vals)
        dlp_pow = np.array(dlp_dat.p_norm_vals)
        dlp_phs = np.array(dlp_dat.phase_rad_vals)
        dlp_frq = np.array(dlp_dat.f_sky_hz_vals)
        geo_spm = np.array(geo_dat.t_oet_spm_vals)

        geo_rho = geo_dat.rho_km_vals
        n_rho = np.size(geo_rho)
        drho = np.zeros(n_rho-1)
        geo_D = geo_dat.D_km_vals
        geo_B = geo_dat.B_deg_vals
        geo_drho = geo_dat.rho_dot_kms_vals
        geo_phi = geo_dat.phi_ora_deg_vals
        t_dlp1 = np.min(dlp_spm)
        t_dlp2 = np.max(dlp_spm)
        t_geo1 = np.min(geo_spm)
        t_geo2 = np.max(geo_spm)

        t1 = np.max([t_dlp1,t_geo1])
        t2 = np.min([t_dlp2,t_geo2])
        if (t1 > t2):
            raise ValueError(
                "Geo and DLP data never overlap. No data available."
            )

        start = np.min((geo_spm >= t1).nonzero())
        finish = np.max((geo_spm <= t2).nonzero())
        tstart = np.min((dlp_spm >= t1).nonzero())
        tfinish = np.max((dlp_spm <= t2).nonzero())
        t_dlp = dlp_spm[tstart:tfinish+1]
        
        for i in range(n_rho-1):
            drho[i] = geo_rho[i+1]-geo_rho[i]
        if not occ:
            if (np.min(drho) < 0.0) and (np.max(drho) > 0.0):
                raise ValueError(
                    "\n\tdrho is positive and negative.\n\
                     \tSet occ to ingress or egress"
                )
            elif (drho > 0).all():
                crange = np.arange(start,finish+1)
            elif (drho < 0).all():
                crange = np.arange(start,finish+1)
                crange = crange[::-1]
            elif (drho == 0).all():
                raise ValueError("drho/dt = 0 for all points.")
            else:
                raise ValueError("drho/dt has invalid values.")
        elif (occ == 'ingress'):
            crange = (drho < 0.0).nonzero()
            if (np.size(crange) == 0):
                raise TypeError("drho is never negative. Use occ = 'egress'")
        elif (occ == 'egress'):
            crange = (drho > 0.0).nonzero()
            if (np.size(crange) == 0):
                raise TypeError("drho is never positive. Use occ = 'ingress'")
        else:
            raise ValueError("Invalid occ keyword: %s" % occ)

        dlp_rho_interp = interpolate.interp1d(geo_spm,geo_rho,kind="linear")
        dlp_rho = dlp_rho_interp(t_dlp)

        rho_min = np.min(dlp_rho)
        rho_max = np.max(dlp_rho)
        self.rho_km_vals = np.arange(rho_min,rho_max,dx)

        drho_interp = interpolate.interp1d(geo_rho,geo_drho,kind="linear")
        self.rho_dot_kms_vals = drho_interp(self.rho_km_vals)
        
        geo_D_interp = interpolate.interp1d(geo_rho,geo_D, kind="linear")
        self.D_km_vals = geo_D_interp(self.rho_km_vals)

        geo_B_interp = interpolate.interp1d(geo_rho, geo_B, kind="linear")
        self.B_rad_vals = np.deg2rad(geo_B_interp(self.rho_km_vals))

        geo_phi_interp = interpolate.interp1d(geo_rho, geo_phi, kind="linear")
        self.phi_rad_vals = np.deg2rad(geo_phi_interp(self.rho_km_vals))

        power_interp = interpolate.interp1d(dlp_rho, dlp_pow, kind="linear")
        self.p_norm_vals = power_interp(self.rho_km_vals)

        phase_interp = interpolate.interp1d(dlp_rho, dlp_phs[tstart:tfinish+1],
                                            kind="linear")
        self.phase_rad_vals = phase_interp(self.rho_km_vals)

        freq_interp = interpolate.interp1d(dlp_rho, dlp_frq[tstart:tfinish+1],
                                           kind="linear")
        self.f_sky_hz_vals = freq_interp(self.rho_km_vals)

        n = np.size(self.rho_km_vals)

        self.t_oet_spm_vals = np.zeros(n)
        self.t_ret_spm_vals = np.zeros(n)
        self.t_set_spm_vals = np.zeros(n)
        self.rho_corr_pole_km_vals = np.zeros(n)
        self.rho_corr_timing_km_vals = np.zeros(n)
        self.phi_rl_rad_vals = np.zeros(n)
        self.raw_tau_threshold_vals = np.zeros(n)

        if verbose:
            print("\tData Extraction Complete.")
        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data": geodata,
            "DLP Data": dlpdata
            }
        input_kwds = {
            "occ":             occ,
            "Use of Verbose":  verbose
            }
        self.history = write_history_dict(input_vars, input_kwds, __file__)

        if verbose:
            print("\tHistory Complete.")
        if verbose:
            print("\tExtract CSV Data Complete.")


class PureCSVReader(object):
    def __init__(self, dat):
        if (not isinstance(dat, str)):
            raise TypeError("Text file must be a string.")
        df = pd.read_csv(dat)
        self.rho_km_vals      = np.array(df.rho_km_vals)
        self.phase_rad_vals   = np.array(df.phase_rad_vals)
        self.p_norm_vals      = np.array(df.p_norm_vals)
        self.phi_rad_vals     = np.array(df.phi_rad_vals)
        self.B_rad_vals       = np.array(df.B_rad_vals)
        self.f_sky_hz_vals    = np.array(df.f_sky_hz_vals)
        self.D_km_vals        = np.array(df.D_km_vals)
        self.rho_dot_kms_vals = np.array(df.rho_dot_kms_vals)

        n_rho = np.size(self.rho_km_vals)
        
        zeros = np.zeros(n_rho)

        # Set fake variables so DiffractionCorrection won't crash.
        self.t_oet_spm_vals = zeros
        self.t_ret_spm_vals = zeros
        self.t_set_spm_vals = zeros
        self.rho_corr_pole_km_vals = zeros
        self.rho_corr_timing_km_vals = zeros
        self.phi_rl_rad_vals = zeros
        self.raw_tau_threshold_vals = zeros
        self.tau_rho = zeros
        self.tau_vals = zeros
        self.phase_vals = zeros
        self.power_vals = zeros

        # Create a history.
        input_vars = {
            "CSV File": dat,
        }

        input_kwds = {}

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        self.rev_info = {
            "rsr_file": "Unknown",
            "band": "Unknown",
            "year": "Unknown",
            "doy": "Unknown",
            "dsn": "Unknown",
            "occ_dir": "Unknown",
            "planetary_occ_flag": "Unknown",
            "rev_num": "Unknown",
            "prof_dir": "Unknown"
        }
