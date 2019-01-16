import numpy as np
import pandas as pd
from scipy import interpolate
from .history import write_history_dict
RADS_PER_DEGS = 0.0174532925199432957692369

def get_geo(geo, verbose=True):
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
    if (not isinstance(geo, str)):
        raise TypeError(
            "\n\tgeo must be a string: '/path/to/geo'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(geo).__name__)
        )
    elif not isinstance(verbose, bool):
        raise TypeError(
            "\n\tverbose must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet verbose=True or verbose=False\n"
            % (type(verbose).__name__)
        )
    else:
        pass

    if verbose:
        print("\tExtracting Geo Data...")

    try:
        dfg = pd.read_csv(geo, delimiter=',',
            names=[
                "t_oet_spm_vals",
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
                "obs_spacecract_lat_deg_vals"
            ]
            )
    except FileNotFoundError:
        raise FileNotFoundError(
            "\n\tYour input geo file does not exists.\n"
            "\tYour file: '%s'\n"
            % (geo)
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
    if (not isinstance(cal, str)):
        raise TypeError(
            "\n\tcal must be a string: '/path/to/cal'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(cal).__name__)
        )
    elif not isinstance(verbose, bool):
        raise TypeError(
            "\n\tverbose must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet verbose=True or verbose=False\n"
            % (type(verbose).__name__)
        )
    else:
        pass

    if verbose:
        print("\tExtracting Cal Data...")

    try:
        dfc = pd.read_csv(cal, delimiter=',',
            names=[
                "spm_vals",
                "f_sky_pred_vals",
                "f_sky_resid_fit_vals",
                "p_free_vals"
                ]
            )
    except FileNotFoundError:
        raise FileNotFoundError(
            "\n\tYour input cal file does not exists.\n"
            "\tYour file: '%s'\n"
            % (cal)
        )

    if verbose:
        print("\tCal Data Complete.")

    return dfc

def get_dlp(dlp, verbose=True):
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
    if (not isinstance(dlp, str)):
        raise TypeError(
            "\n\tdlp must be a string: '/path/to/dlp'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(dlp).__name__)
        )
    elif not isinstance(verbose, bool):
        raise TypeError(
            "\n\tverbose must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet verbose=True or verbose=False\n"
            % (type(verbose).__name__)
        )
    else:
        pass

    if verbose:
        print("\tExtracting DLP Data...")

    try:
        dfd = pd.read_csv(
            dlp, delimiter=',',
            names=[
                "rho_km_vals",
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
                "B_deg_vals"
            ]
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            "\n\tYour input dlp file does not exists.\n"
            "\tYour file: '%s'\n"
            % (dlp)
        )

    if verbose:
        print("\tDLP Data Complete")
    return dfd

def get_tau(tau, verbose=True):
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
    if (not isinstance(tau, str)):
        raise TypeError(
            "\n\ttau must be a string: '/path/to/tau'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(tau).__name__)
        )
    elif not isinstance(verbose, bool):
        raise TypeError(
            "\n\tverbose must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet verbose=True or verbose=False\n"
            % (type(verbose).__name__)
        )
    else:
        pass

    if verbose:
        print("\tExtracting Tau Data...")

    try:
        dft = pd.read_csv(tau, delimiter=',',
            names=[
                "rho_km_vals",
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
                "B_deg_vals"
            ]
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            "\n\tYour input tau file does not exists.\n"
            "\tYour file: '%s'\n"
            % (tau)
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
        Dependencies:
            #. pandas
            #. numpy
            #. scipy
            #. rss_ringoccs
    """
    def __init__(self, geo, cal, dlp, tau=None, verbose=True):
        if (not isinstance(geo, str)):
            raise TypeError(
                "geo must be a string: '/path/to/geo'\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: str\n"
                % (type(geo).__name__)
            )
        elif (not isinstance(cal, str)):
            raise TypeError(
                "cal must be a string: '/path/to/cal'\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: str\n"
                % (type(cal).__name__)
            )
        elif (not isinstance(dlp, str)):
            raise TypeError(
                "dlp must be a string: '/path/to/dlp'\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: str\n"
                % (type(dlp).__name__)
            )
        elif not isinstance(verbose, bool):
            raise TypeError(
                "\n\tverbose must be Boolean: True/False\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: bool\n"
                "\tSet verbose=True or verbose=False\n"
                % (type(verbose).__name__)
            )
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
        geo_dat = get_geo(self.geo, verbose=verbose)
        cal_dat = get_cal(self.cal, verbose=verbose)
        dlp_dat = get_dlp(self.dlp, verbose=verbose)

        if verbose:
            print("\tRetrieving Variables...")

        try:
            # Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"

            # Ring radius.
            self.rho_km_vals = np.array(dlp_dat.rho_km_vals)
            if not (np.all(np.isreal(self.rho_km_vals))):
                raise ValueError(
                    "\n\t\trho_km_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.rho_km_vals = self.rho_km_vals.astype(float)

            if (np.min(self.rho_km_vals) <= 0.0):
                raise ValueError(
                    "\n\t\trho_km_vals has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tdlp CSV for errors.\n"
                    "\t\t\tMinimum rho value: %f"
                    % np.min(self.rho_km_vals)
                )
            else:
                pass

            # Raw optical depth.
            errmess = "raw_tau_vals"
            raw_tau_vals = np.array(dlp_dat.raw_tau_vals)
            if not (np.all(np.isreal(raw_tau_vals))):
                raise ValueError(
                    "\n\t\traw_tau_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                raw_tau_vals = raw_tau_vals.astype(float)

            # Phase of signal.
            errmess = "phase_deg_vals"
            phase_deg_vals = np.array(dlp_dat.phase_deg_vals)
            if not (np.all(np.isreal(phase_deg_vals))):
                raise ValueError(
                    "\n\t\tphase_deg_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                phase_deg_vals = phase_deg_vals.astype(float)

            if (np.max(np.abs(phase_deg_vals)) >= 360+1.e-4):
                raise ValueError(
                    "\n\t\tphase_deg_vals has values beyond\n"
                    "\t\t360 degrees. Please check your dlp\n"
                    "\t\tCSV for errors.\n"
                    "\t\t\tMaximum of |phase|: %f"
                    % (np.max(np.abs(phase_deg_vals)))
                )
            else:
                pass

            # Ring azimuth angle.
            errmess = "phi_ora_deg_vals"
            phi_ora_deg_vals = np.array(geo_dat.phi_ora_deg_vals)
            if not (np.all(np.isreal(phi_ora_deg_vals))):
                raise ValueError(
                    "\n\t\tphi_ora_deg_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                phi_ora_deg_vals = phi_ora_deg_vals.astype(float)

            if (np.max(np.abs(phi_ora_deg_vals)) >= 360+1.e-4):
                raise ValueError(
                    "\n\t\tphi_ora_deg_vals has values beyond\n"
                    "\t\t360 degrees. Please check your geo\n"
                    "\t\tCSV for errors.\n"
                    "\t\t\tMaximum of |phi|: %f"
                    % (np.max(np.abs(phi_ora_deg_vals)))
                )
            else:
                pass

            # Ring opening angle.
            errmess = "B_deg_vals"
            B_deg_vals = np.array(geo_dat.B_deg_vals)
            if not (np.all(np.isreal(B_deg_vals))):
                raise ValueError(
                    "\n\t\tB_deg_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                B_deg_vals = B_deg_vals.astype(float)

            if (np.max(np.abs(B_deg_vals)) >= 360+1.e-4):
                raise ValueError(
                    "\n\t\tB_deg_vals has values beyond\n"
                    "\t\t360 degrees. Please check your geo\n"
                    "\t\tCSV for errors.\n"
                    "\t\t\tMaximum of |B|: %f"
                    % (np.max(np.abs(B_deg_vals)))
                )
            else:
                pass

            # Retrieve time variables (Earth, Ring, and Spacecraft ET).
            ermess = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(dlp_dat.t_oet_spm_vals)
            if not (np.all(np.isreal(self.t_oet_spm_vals))):
                raise ValueError(
                    "\n\t\tt_oet_spm_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.t_oet_spm_vals = self.t_oet_spm_vals.astype(float)

            ermess = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(dlp_dat.t_ret_spm_vals)
            if not (np.all(np.isreal(self.t_ret_spm_vals))):
                raise ValueError(
                    "\n\t\tt_ret_spm_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.t_ret_spm_vals = self.t_ret_spm_vals.astype(float)

            ermess = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(dlp_dat.t_set_spm_vals)
            if not (np.all(np.isreal(self.t_set_spm_vals))):
                raise ValueError(
                    "\n\t\tt_set_spm_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.t_set_spm_vals = self.t_set_spm_vals.astype(float)

            # Pole correction for rho.
            ermess = "rho_corr_pole_km_vals"
            self.rho_corr_pole_km_vals = np.array(dlp_dat.rho_corr_pole_km_vals)
            if not (np.all(np.isreal(self.rho_corr_pole_km_vals))):
                raise ValueError(
                    "\n\t\trho_corr_pole_km_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.rho_corr_pole_km_vals = self.rho_corr_pole_km_vals.astype(float)

            # Timing corrections in ring radius.
            ermess = "rho_corr_pole_km_vals"
            self.rho_corr_timing_km_vals = np.array(dlp_dat.rho_corr_timing_km_vals)
            if not (np.all(np.isreal(self.rho_corr_timing_km_vals))):
                raise ValueError(
                    "\n\t\trho_corr_timing_km_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals.astype(float)

            # Ring longitude angle.
            errmess = "phi_rl_deg_vals"
            self.phi_rl_deg_vals = np.array(geo_dat.phi_rl_deg_vals)
            if not (np.all(np.isreal(self.phi_rl_deg_vals))):
                raise ValueError(
                    "\n\t\tphi_rl_deg_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                self.phi_rl_deg_vals = self.phi_rl_deg_vals.astype(float)

            if (np.max(np.abs(self.phi_rl_deg_vals)) >= 360+1.e-4):
                raise ValueError(
                    "\n\t\tphi_rl_deg_vals has values beyond\n"
                    "\t\t360 degrees. Please check your geo\n"
                    "\t\tCSV for errors.\n"
                    "\t\t\tMaximum of |phi|: %f"
                    % (np.max(np.abs(self.phi_rl_deg_vals)))
                )
            else:
                pass

            # Optical depth threshold of diffraction profile.
            errmess = "raw_tau_threshold_vals"
            self.raw_tau_threshold_vals = np.array(dlp_dat.raw_tau_threshold_vals)
            if not (np.all(np.isreal(self.raw_tau_threshold_vals))):
                raise ValueError(
                    "\n\t\traw_tau_threshold_vals is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your dlp CSV for errors.\n"
                )
            else:
                self.raw_tau_threshold_vals = self.raw_tau_threshold_vals.astype(float)

            # Ring radius from GEO file.
            errmess = "geo_rho"
            geo_rho = np.array(geo_dat.rho_km_vals)
            if not (np.all(np.isreal(geo_rho))):
                raise ValueError(
                    "\n\t\tgeo_rho is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                geo_rho = geo_rho.astype(float)

            if (np.min(geo_rho) <= 0.0):
                raise ValueError(
                    "\n\t\tgeo_rho has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tgeo CSV for errors.\n"
                    "\t\t\tMinimum rho value: %f"
                    % np.min(geo_rho)
                )
            else:
                pass

            # RIP-Spacecraft disstance from GEO file.
            errmess = "geo_D"
            geo_D = np.array(geo_dat.D_km_vals)
            if not (np.all(np.isreal(geo_D))):
                raise ValueError(
                    "\n\t\tgeo_D is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                geo_D = geo_D.astype(float)

            if (np.min(geo_D) <= 0.0):
                raise ValueError(
                    "\n\t\tgeo_D has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tgeo CSV for errors.\n"
                    "\t\t\tMinimum D value: %f"
                    % np.min(geo_D)
                )
            else:
                pass

            # drho/dt from GEO file.
            errmess = "geo_drho"
            geo_drho = np.array(geo_dat.rho_dot_kms_vals)
            if not (np.all(np.isreal(geo_D))):
                raise ValueError(
                    "\n\t\tgeo_drho is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your geo CSV for errors.\n"
                )
            else:
                geo_drho = geo_drho.astype(float)

            # Sky frequency
            errmess = "f_sky_pred"
            f_sky_pred = np.array(cal_dat.f_sky_pred_vals)
            if not (np.all(np.isreal(f_sky_pred))):
                raise ValueError(
                    "\n\t\tf_sky_pred is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your cal CSV for errors.\n"
                )
            else:
                f_sky_pred = f_sky_pred.astype(float)

            if (np.min(f_sky_pred) <= 0.0):
                raise ValueError(
                    "\n\t\tf_sky_pred has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tcal CSV for errors.\n"
                    "\t\t\tMinimum frequency value: %f"
                    % np.min(f_sky_pred)
                )
            else:
                pass

            errmess = "f_sky_resid"
            f_sky_resid = np.array(cal_dat.f_sky_resid_fit_vals)
            if not (np.all(np.isreal(f_sky_resid))):
                raise ValueError(
                    "\n\t\tf_sky_resid is not an array of real\n"
                    "\t\tvalued floating point numbers. Please\n"
                    "\t\tcheck your cal CSV for errors.\n"
                )
            else:
                f_sky_resid = f_sky_resid.astype(float)

            if (np.min(f_sky_resid) <= 0.0):
                raise ValueError(
                    "\n\t\tf_sky_resid has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tcal CSV for errors.\n"
                    "\t\t\tMinimum frequency value: %f"
                    % np.min(f_sky_resid)
                )
            else:
                pass

            errmess = "f_sky_resid"
            f_sky_raw_vals = f_sky_pred - f_sky_resid
            if (np.min(f_sky_raw_vals) <= 0.0):
                raise ValueError(
                    "\n\t\tf_sky_raw_vals has values that are\n"
                    "\t\tnot positive. Please check your input\n"
                    "\t\tcal CSV for errors.\n"
                    "\t\t\tMinimum frequency value: %f"
                    % np.min(f_sky_raw_vals)
                )
            else:
                pass
        except (ValueError, TypeError, NameError, AttributeError) as err:
            raise TypeError(
                "\n\tError occured while trying to extract\n"
                "\tthe variable: %s\n"
                "\tOriginal error message:\n"
                "%s"
                % (errmess, err)
            )

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
                "\n\tdrho/dt has positive and negative values.\n"
                "\tPlease check your dlp CSV for errors."
            )
        elif (np.size((drdt == 0).nonzero()) != 0):
            raise ValueError(
                "\n\tdrho/dt has zero valued elements.\n"
                "\tPlease check your dlp CSV for errors."
            )
        elif (drdt < 0.0).all():
            occ = 'ingress'
        elif (drdt > 0.0).all():
            occ = 'egress'
        else:
            raise ValueError("\n\tBad DLP: drho/dt has incompatible elements")

        if (occ == 'ingress'):
            crange = (geo_drho < 0.0).nonzero()
        elif (occ == 'egress'):
            crange = (geo_drho > 0.0).nonzero()
        else:
            crange_e = (geo_drho > 0.0).nonzero()
            crange_i = (geo_drho < 0.0).nonzero()
            n_e = np.size(crange_e)
            n_i = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    "\n\trho_dot_kms_vals has positive and negative values.\n"
                    "\tThis is likely a chord occultation. Set occ='ingress'\n"
                    "\tto examine the ingress portion, and occ='egress'\n"
                    "\tfor the egress porition."
                )
            elif (n_e == 0) and (n_i == 0):
                raise ValueError(
                    "\n\trho_dot_kms_vals is either empty or zero.\n"
                    "\tPlease check your geo dlp CSV for errors."
                )
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else:
                raise TypeError(
                    "\n\tCould not determine what type of\n"
                    "\toccultation this is. Please check you\n"
                    "geo CSV for errors."
                )

            del n_e, n_i, crange_e, crange_i

        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else:
                raise ValueError("Bad occ input: Set 'egress' or 'ingress'")

            raise ValueError(
                "\n\tBad occ Input.\n"
                "\t%s\n"
                "\tOccultation Type: %s" % (mes, occ)
            )

        if verbose:
            print("\tInterpolating Data...")

        geo_drho = geo_drho[crange]
        geo_rho = geo_rho[crange]
        geo_D = geo_D[crange]
        phi_ora_deg_vals = phi_ora_deg_vals[crange]
        B_deg_vals = B_deg_vals[crange]
        rmin = np.min(geo_rho)
        rmax = np.max(geo_rho)
        rfin = int(np.max((rmax-self.rho_km_vals>=0.0).nonzero()))
        rstart = int(np.min((self.rho_km_vals-rmin>=0.0).nonzero()))
        n_rho_vals = np.size(self.rho_km_vals)
        n_f_vals = np.size(f_sky_raw_vals)
        frange = np.arange(n_f_vals)
        xrange = np.arange(n_rho_vals)*(n_f_vals-1.0)/(n_rho_vals-1.0)
        interp = interpolate.interp1d(geo_rho, geo_D, kind='cubic')
        self.D_km_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, geo_drho, kind='cubic')
        self.rho_dot_kms_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, phi_ora_deg_vals, kind="cubic")
        self.phi_rad_vals = interp(self.rho_km_vals)*RADS_PER_DEGS
        interp = interpolate.interp1d(geo_rho, B_deg_vals, kind="cubic")
        self.B_rad_vals = interp(self.rho_km_vals)*RADS_PER_DEGS
        interp = interpolate.interp1d(geo_rho, phi_ora_deg_vals, kind="cubic")
        self.phi_rl_rad_vals = interp(self.rho_km_vals)*RADS_PER_DEGS
        interp = interpolate.interp1d(frange, f_sky_raw_vals, kind='cubic')
        self.f_sky_hz_vals = interp(xrange)
        raw_mu = np.sin(np.abs(self.B_rad_vals))
        self.p_norm_vals = np.exp(-raw_tau_vals/raw_mu)
        self.t_ret_spm_vals = self.t_ret_spm_vals[rstart:rfin+1]
        self.t_set_spm_vals = self.t_set_spm_vals[rstart:rfin+1]
        self.t_oet_spm_vals = self.t_oet_spm_vals[rstart:rfin+1]
        self.phase_rad_vals = phase_deg_vals[rstart:rfin+1]*RADS_PER_DEGS
        self.rho_km_vals = self.rho_km_vals[rstart:rfin+1]
        self.p_norm_vals = self.p_norm_vals[rstart:rfin+1]
        self.rho_corr_pole_km_vals = self.rho_corr_pole_km_vals[rstart:rfin+1]
        self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals[rstart:rfin+1]

        del f_sky_raw_vals, rmin, rmax, rstart, rfin, n_rho_vals, n_f_vals
        del interp, frange, xrange, phi_ora_deg_vals, raw_tau_vals
        del phase_deg_vals, raw_mu, B_deg_vals, geo_rho, geo_D
        del geo_drho, crange, geo_dat, cal_dat, dlp_dat

        if (not isinstance(tau, type(None))):
            if (not isinstance(tau, str)):
                raise TypeError("taudata must be a string: '/path/to/taudata'")
            else:
                tau_dat = get_tau(self.tau, verbose=verbose)
                tp = tau_dat.phase_deg_vals*RADS_PER_DEGS
                tm = tau_dat.B_deg_vals*RADS_PER_DEGS
                tm = np.sin(np.abs(tm))
                tr = tau_dat.rho_km_vals
                tt = tau_dat.raw_tau_vals
                rmin = np.min(tr)
                rmax = np.max(tr)
                rfin = int(np.max((rmax-self.rho_km_vals>=0).nonzero()))
                rstart = int(np.min((self.rho_km_vals-rmin>=0).nonzero()))
                self.tau_rho = self.rho_km_vals[rstart:rfin+1]
                interp = interpolate.interp1d(tr, tt, kind='cubic')
                self.tau_vals = interp(self.tau_rho)
                interp = interpolate.interp1d(tr, tp, kind='cubic')
                self.phase_vals = interp(self.tau_rho)
                interp = interpolate.interp1d(tr, tm, kind='cubic')
                tm = interp(self.tau_rho)
                self.power_vals = np.exp(-self.tau_vals/tm)
                del tau_dat, tp, tm, tr, tt, rmin, rmax, rfin, rstart, interp
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
        try:
            var = var.split("_")
            band = var[3][0]
            year = var[1]
            doy = var[2]
            dsn = "DSS-%s" % (var[3][1:])
            occ_dir = var[4]
            rev_num = "Unknown"
            if (occ_dir == "E"):
                prof_dir = "EGRESS"
            else:
                occ_dir = "INGRESS"
        except:
            var = "Unknown"
            band = "Unknown"
            year = "Unknown"
            doy = "Unknown"
            dsn = "Unknown"
            occ_dir = "Unknown"
            rev_num = "Unknown"
            prof_dir = "Unknown"

        self.rev_info = {
            "rsr_file": "Unknown",
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
