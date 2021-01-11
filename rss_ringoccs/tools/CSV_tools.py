"""
    :Purpose:
            Provide tools for reading in .TAB and .CSV files and converting
            the data into a usable instance of the DLP class.

    :Dependencies:
        #. pandas
        #. numpy
"""

import numpy as np
import pandas as pd
from .history import write_history_dict, date_to_rev, rev_to_occ_info
from . import error_check

def get_geo(geo, verbose=True, use_deprecate=False):
    """
    To extract a pandas DataFrame from a given GEO.TAB or GEO.CSV file.

    Arguments
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

    Keywords
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
                                     "vz_kms_vals"],
                                     usecols=list(range(18)))
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
                                     "obs_spacecract_lat_deg_vals"],
                                     usecols=list(range(19)))
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
                                 "p_free_vals"],
                                 usecols=list(range(4)))
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
                                     "B_deg_vals"],
                                     usecols=list(range(12)))
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
                                     "B_deg_vals"], usecols=list(range(13)))
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
            # Read in DLP. Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = np.array(dlp_dat.rho_km_vals)
            errmess = "raw_tau_vals"
            self.raw_tau_vals = np.array(dlp_dat.raw_tau_vals)
            errmess = "phase_rad_vals"
            self.phase_rad_vals = np.radians(np.array(dlp_dat.phase_deg_vals))
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = np.radians(np.array(dlp_dat.phi_ora_deg_vals))
            errmess = "B_rad_vals"
            self.B_rad_vals = np.radians(np.array(dlp_dat.B_deg_vals))
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
            self.phi_rl_rad_vals = np.radians(np.array(dlp_dat.phi_rl_deg_vals))
            errmess = "raw_tau_threshold_vals"
            self.raw_tau_threshold_vals = np.array(dlp_dat.raw_tau_threshold_vals)

            # Grab info from CAL
            errmess = "f_sky_pred"
            f_sky_pred = np.array(cal_dat.f_sky_pred_vals)
            errmess = "f_sky_resid"
            f_sky_resid = np.array(cal_dat.f_sky_resid_fit_vals)
            errmess = "f_sky_resid"
            self.f_sky_hz_vals = f_sky_pred - f_sky_resid

            # Grab remaining necessary info from GEO
            errmess = "geo_rho"
            geo_rho = np.array(geo_dat.rho_km_vals)
            errmess = "D_km_vals"
            self.D_km_vals = np.array(geo_dat.D_km_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(geo_dat.rho_dot_kms_vals)
            rx_geo = np.array(geo_dat.rx_km_vals)
            ry_geo = np.array(geo_dat.ry_km_vals)
            rz_geo = np.array(geo_dat.rz_km_vals)

            errmess = "cal_spm"
            cal_spm = np.array(cal_dat.spm_vals)
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

        error_check.check_lengths(self.t_set_spm_vals, self.rho_km_vals,
                                  "t_set_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.t_ret_spm_vals, self.rho_km_vals,
                                  "t_ret_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.t_oet_spm_vals, self.rho_km_vals,
                                  "t_oet_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phase_rad_vals, self.rho_km_vals,
                                  "phase_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phi_rad_vals, self.rho_km_vals,
                                  "phi_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.raw_tau_vals, self.rho_km_vals,
                                  "raw_tau_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.B_rad_vals, self.rho_km_vals,
                                  "B_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.D_km_vals, geo_rho,
                                  "D_km_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.rho_dot_kms_vals, geo_rho,
                                  "D_km_vals", "rho_km_vals", fname)

        if verbose:
            print("\tComputing Variables...")

        dr = np.diff(self.rho_km_vals)
        dt = np.diff(self.t_set_spm_vals)

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

            # The rev is ingress, so flip GEO variables.
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])
            self.D_km_vals = self.D_km_vals[::-1]
            geo_rho = geo_rho[::-1]
        elif (drdt > 0.0).all():
            pass
        else:
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not determine occultation type.
                """ % (fname)
            )

        if verbose:
            print("\tInterpolating Data...")

        del geo_dat, cal_dat, dlp_dat

        raw_mu = np.sin(np.abs(self.B_rad_vals))
        self.p_norm_vals = np.exp(-self.raw_tau_vals/raw_mu)

        self.f_sky_hz_vals = np.interp(self.t_oet_spm_vals, cal_spm,
                                       self.f_sky_hz_vals)
        self.D_km_vals = np.interp(self.rho_km_vals, geo_rho, self.D_km_vals)
        self.rho_dot_kms_vals = np.interp(self.rho_km_vals, geo_rho,
                                          self.rho_dot_kms_vals)
        self.rx_km_vals = np.interp(self.rho_km_vals, geo_rho, rx_geo)
        self.ry_km_vals = np.interp(self.rho_km_vals, geo_rho, ry_geo)
        self.rz_km_vals = np.interp(self.rho_km_vals, geo_rho, rz_geo)

        del raw_mu, geo_rho

        if (not isinstance(tau, type(None))):
            tau_dat = get_tau(self.tau, verbose=verbose, use_deprecate=use_deprecate)
            tm = np.sin(np.abs(np.deg2rad(tau_dat.B_deg_vals)))
            rmin = np.min(tau_dat.rho_km_vals)
            rmax = np.max(tau_dat.rho_km_vals)
            rfin = int(np.max((rmax-self.rho_km_vals>=0).nonzero()))
            rstart = int(np.min((self.rho_km_vals-rmin>=0).nonzero()))
            self.tau_rho = self.rho_km_vals[rstart:rfin+1]
            self.tau_vals = np.interp(self.tau_rho, tau_dat.rho_km_vals,
                                      tau_dat.raw_tau_vals)
            self.phase_vals = np.interp(self.tau_rho, tau_dat.rho_km_vals,
                                        np.deg2rad(tau_dat.phase_deg_vals))
            tm = np.interp(self.tau_rho, tau_dat.rho_km_vals, tm)
            self.power_vals = np.exp(-self.tau_vals/tm)
            del tau_dat, tm, rmin, rmax, rfin, rstart
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

        try:
            var = (geo.split("/")[-1]).split("_")
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
        except:
            print(
                """
                    \r\tError: rss_ringoccs
                    \r\t\ttools.CSV_tools.ExtractCSVData\n
                    \r\tCould not set rev_info. Returning data without this.
                    \r\twrite_file option will not work with this instance of
                    \r\tthe ExtractCSVData class. This can still be used with
                    \r\tDiffractionCorrection if write_file=False is set.\n
                    \r\tTo correct, re-run within pipeline directory.\n
                    \rOriginal error message:\n%s
                """ % err
            )

            self.rev_info = None

        if verbose:
            print("\tHistory Complete.")
            print("\tExtract CSV Data Complete.")


class GetUranusData(object):
    def __init__(self, geodata, dlpdata, tau=None, verbose=False):
        fname = "tools.CSV_tools.ExtractCSVData"
        error_check.check_type(geodata, str, "geo", fname)
        error_check.check_type(dlpdata, str, "dlp", fname)
        error_check.check_type(verbose, bool, "verbose", fname)
        if (not isinstance(tau, type(None))):
            error_check.check_type(tau, str, "tau", fname)

        geo_dat = get_geo(geodata, verbose=verbose)

        dlp_dat = pd.read_csv(
            dlpdata, delimiter=',',
            names=[
                "rho_km_vals",
                "null1",
                "null2",
                "phi_rl_rad_vals",
                "phi_rad_vals",
                "p_norm_vals",
                "raw_tau_vals",
                "phase_rad_vals",
                "raw_tau_threshold_vals",
                "t_oet_spm_vals",
                "t_ret_spm_vals",
                "t_set_spm_vals",
                "B_rad_vals",
                "rho_dot_kms_vals",
                "F_km_vals",
                "D_km_vals",
                "f_sky_hz_vals"
            ],
            usecols=list(range(17))
        )

        if verbose:
            print("\tRetrieving Variables...")

        try:
            # Read in DLP. Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = np.array(dlp_dat.rho_km_vals)
            errmess = "raw_tau_vals"
            self.raw_tau_vals = np.array(dlp_dat.raw_tau_vals)
            errmess = "phase_rad_vals"
            self.phase_rad_vals = np.array(dlp_dat.phase_rad_vals)
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = np.array(dlp_dat.phi_rad_vals)
            errmess = "B_rad_vals"
            self.B_rad_vals = np.array(dlp_dat.B_rad_vals)
            errmess = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(dlp_dat.t_oet_spm_vals)
            errmess = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(dlp_dat.t_ret_spm_vals)
            errmess = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(dlp_dat.t_set_spm_vals)
            errmess = "rho_corr_pole_km_vals"
            self.rho_corr_pole_km_vals = np.array(dlp_dat.rho_km_vals)
            errmess = "rho_corr_pole_km_vals"
            self.rho_corr_timing_km_vals = np.array(dlp_dat.rho_km_vals)
            errmess = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = np.array(dlp_dat.phi_rl_rad_vals)
            errmess = "raw_tau_threshold_vals"
            self.raw_tau_threshold_vals = np.array(dlp_dat.raw_tau_threshold_vals)
            errmess = "p_norm_vals"
            self.p_norm_vals = np.array(dlp_dat.p_norm_vals)
            errmess = "D_km_vals"
            self.D_km_vals = np.array(dlp_dat.D_km_vals)
            errmess = "f_sky_hz_vals"
            self.f_sky_hz_vals = np.array(dlp_dat.f_sky_hz_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(dlp_dat.rho_dot_kms_vals)

            # Grab remaining necessary info from GEO
            errmess = "geo_rho"
            geo_rho = np.array(geo_dat.rho_km_vals)
            errmess = "rx_geo"
            rx_geo = np.array(geo_dat.rx_km_vals)
            errmess = "ry_geo"
            ry_geo = np.array(geo_dat.ry_km_vals)
            errmess = "rz_geo"
            rz_geo = np.array(geo_dat.rz_km_vals)

        except (ValueError, TypeError, NameError, AttributeError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not convert %s into a numpy array.
                    \r\tCheck input GEO, CAL, and DLP files.
                """ % (fname, errmess)
            )

        error_check.check_is_real(self.raw_tau_threshold_vals,
                                  "raw_tau_threshold_vals", fname)
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

        error_check.check_lengths(self.t_set_spm_vals, self.rho_km_vals,
                                  "t_set_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.t_ret_spm_vals, self.rho_km_vals,
                                  "t_ret_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.t_oet_spm_vals, self.rho_km_vals,
                                  "t_oet_spm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phase_rad_vals, self.rho_km_vals,
                                  "phase_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phi_rad_vals, self.rho_km_vals,
                                  "phi_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.raw_tau_vals, self.rho_km_vals,
                                  "raw_tau_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.B_rad_vals, self.rho_km_vals,
                                  "B_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.D_km_vals, self.rho_km_vals,
                                  "D_km_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.rho_dot_kms_vals, self.rho_km_vals,
                                  "D_km_vals", "rho_km_vals", fname)

        if verbose:
            print("\tComputing Variables...")

        dr = np.diff(self.rho_km_vals)
        dt = np.diff(self.t_set_spm_vals)

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

            # The rev is ingress, so flip GEO variables.
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])
            self.D_km_vals = self.D_km_vals[::-1]
            geo_rho = geo_rho[::-1]
        elif (drdt > 0.0).all():
            pass
        else:
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not determine occultation type.
                """ % (fname)
            )

        if verbose:
            print("\tInterpolating Data...")

        self.rx_km_vals = np.interp(self.rho_km_vals, geo_rho, rx_geo)
        self.ry_km_vals = np.interp(self.rho_km_vals, geo_rho, ry_geo)
        self.rz_km_vals = np.interp(self.rho_km_vals, geo_rho, rz_geo)

        if (not isinstance(tau, type(None))):
            tau_dat = get_tau(tau, verbose=verbose, use_deprecate=False)
            tm = np.sin(np.abs(np.deg2rad(tau_dat.B_deg_vals)))
            rmin = np.min(tau_dat.rho_km_vals)
            rmax = np.max(tau_dat.rho_km_vals)
            rfin = int(np.max((rmax-self.rho_km_vals>=0).nonzero()))
            rstart = int(np.min((self.rho_km_vals-rmin>=0).nonzero()))
            self.tau_rho = self.rho_km_vals[rstart:rfin+1]
            self.tau_vals = np.interp(self.tau_rho, tau_dat.rho_km_vals,
                                      tau_dat.raw_tau_vals)
            self.phase_vals = np.interp(self.tau_rho, tau_dat.rho_km_vals,
                                        np.deg2rad(tau_dat.phase_deg_vals))
            tm = np.interp(self.tau_rho, tau_dat.rho_km_vals, tm)
            self.power_vals = np.exp(-self.tau_vals/tm)
            del tau_dat, tm, rmin, rmax, rfin, rstart
        else:
            self.tau_rho = None
            self.tau_vals = None
            self.phase_vals = None
            self.power_vals = None

        if verbose:
            print("\tData Extraction Complete.")
            print("\tWriting History...")

        input_vars = {
            "GEO Data": geodata,
            "DLP Data": dlpdata
            }
        input_kwds = {
            "Use of Verbose":  verbose
            }
        self.history = write_history_dict(input_vars, input_kwds, __file__)

        if verbose:
            printf("Data Extraction Complete.")

