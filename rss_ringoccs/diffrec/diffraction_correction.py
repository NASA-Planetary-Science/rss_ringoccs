"""
    Purpose:
        Provide the DiffractionCorrection class for
        performing the necessary mathematics to correct
        for the diffraction effects that are obtained
        during occultation observations of planetary
        rings using radio waves.
    Dependencies:
        #. numpy
        #. scipy
        #. rss_ringoccs
"""
# Import dependencies for the diffcorr module
import numpy as np
import sys
from rss_ringoccs.tools.history import write_history_dict
from rss_ringoccs.tools.write_output_files import write_output_files
from rss_ringoccs.tools import error_check
from . import special_functions, window_functions
try:
    from rss_ringoccs._ufuncs import _diffraction_functions
except:
    print(
        """
            Error: rss_ringoccs.diffrec.diffraction_correction
            \tCould Not Import C Code. Stricly Using Python Code.
            \tThis is signicantly slower. There was most likely an error
            \tin your installation of rss_ringoccs. To use the C Code,
            \tdownload a C Compiler (GCC) and see the User's Guide for
            \tinstallation instructions.
        """
    )

# Declare constant for the speed of light (km/s)
SPEED_OF_LIGHT_KM = 299792.4580

# Declare constants for multiples of pi.
HALF_PI = 1.570796326794896619231322
TWO_PI = 6.283185307179586476925287

# Dictionary containing regions of interest within the Saturnian Rings.
region_dict = {
    'all':             [1.0, 400000.0],
    'besselbarnard':   [120210.0, 120330.0],
    'bessel-barnard':  [120210.0, 120330.0],
    'cringripples':    [77690.0, 77760.0],
    'encke':           [132900.0, 134200.0],
    'enckegap':        [132900.0, 134200.0],
    'herschel':        [118100.0, 118380.0],
    'herschelgap':     [118100.0, 118380.0],
    'huygens':         [117650.0, 117950.0],
    'huygensringlet':  [117650.0, 117950.0],
    'janusepimetheus': [96200.0, 96800.0],
    'jeffreys':        [118900.0, 119000.0],
    'jeffreysgap':     [118900.0, 119000.0],
    'kuiper':          [119300.0, 119500.0],
    'kuipergap':       [119300.0, 119500.0],
    'maxwell':         [87410.0, 87610.0],
    'maxwellringlet':  [87410.0, 87610.0],
    'russell':         [118550.0, 118660.0],
    'russellgap':      [118550.0, 118660.0],
    'titan':           [77870.0, 77930.0],
    'titanringlet':    [77870.0, 77930.0]
}


class DiffractionCorrection(object):
    """
        Purpose:
            Perform diffraction correction for a ring occultation
            on a data set that is a near radially symmetric function
            of the ring radius, or ring intercept point (RIP).
        Arguments:
            :DLP (*object*):
                The data set, usually an instance of the
                DiffractionLimitedProfile class from the rss_ringoccs
                Calibration subpackage. This instance MUST contain
                the following attributes and have the same names.

                |   rho_km_vals:      Ring Radius (km)
                |   phi_rad_vals:     Ring Azimuth Angle (Radians)
                |   p_norm_vals:      Normalized Power
                |   phase_rad_vals:   Phase (Radians)
                |   B_rad_vals:       Elevation Angle (Radians)
                |   D_km_vals:        RIP-Distance (km)
                |   f_sky_hz_vals:    Sky Frequency (Hertz)
                |   rho_dot_kms_vals: RIP-velocity (km/s)
                |   history:          History dictionary

            :res (*float* or *int*):
                The requested resolution for processing (km). This
                must be a positive real number.
        Keywords:
            :rng (*list* or *str*):
                The requested range for diffraction correction.
                Preferred input is rng = [a,b]. Arrays are
                allowed and the range will be set as:

                |    rng = [MIN(array), MAX(array)]

                Finally, certain strings containing a few of the
                regions of interests within the rings of Saturn
                are allowed. Permissable strings are:

                |    'all'             [1.0, 400000.0]
                |    'cringripples'    [77690.0, 77760.0]
                |    'encke'           [132900.0, 134200.0]
                |    'enckegap'        [132900.0, 134200.0]
                |    'janusepimetheus' [96200.0, 96800.0]
                |    'maxwell'         [87410.0, 87610.0]
                |    'maxwellringlet'  [87410.0, 87610.0]
                |    'titan'           [77870.0, 77930.0]
                |    'titanringlet'    [77870.0, 77930.0]
                |    'huygens'         [117650.0, 117950.0]
                |    'huygensringlet'  [117650.0, 117950.0]

                Strings are neither case nor space sensitive.
                For other planets use rng = [a,b]. Default value
                is set to 'all' which processes [1, 400000]
                Values MUST be set in kilometers.
            :wtype (*str):
                The requested tapering function for diffraction
                correction. A string with several allowed inputs:

                |    'rect'      Rectangular Window.
                |    'coss'      Squared Cosine Window.
                |    'kb20'      Kaiser-Bessel 2.0 Window.
                |    'kb25'      Kaiser-Bessel 2.5 Window.
                |    'kb35'      Kaiser-Bessel 3.5 Window.
                |    'kbmd20'    Modified kb20 Window.
                |    'kbmd25'    Modified kb25 Window.

                The variable is neither case nor space sensitive.
                Default window is set to 'kb25'. See window_functions
                submodule for further documentation.
            :fwd (*bool*):
                A Boolean for determining whether or not
                forward modelling will be computed. This is good
                starting point for deciding if the diffraction
                correction is physically significant or valid. If
                the reconstruction is good, the forward model
                should reproduce the p_norm_vals attribute from
                the input DLP instance. Default is set to False.
            :norm (*bool*):
                A Boolean for determining whether or not the
                reconstructed complex transmittance is normalize
                by the window width. This normalization is the
                complex transmittance that is computed by using
                free space divided by the complex transmittance
                that is computed using free space weighted by the
                selected tapering function. Default is True.
            :bfac (*bool*):
                A Boolean for determining whether or not the
                'b' factor in the window width computation is
                used. This is equivalent to setting the Allen
                Deviation for the spacecraft to a positive value
                or to zero. If set to False, the Allen Deviation
                is assumed to be zero. If set to True the Allen
                Deviation is set to 2e-13, or whichever number you
                wish to specify in the sigma keyword (See below).
                Default is True.
            :sigma (*float*):
                The Allen deviation for the spacecraft. If the bfac
                keyword (See above) is set to False, this is ignored.
                If bfac is set to True, and sigma is NOT specified,
                then sigma=2e-13 will be used, which is the Allen
                deviation for Cassini with 1 second integration time.
                For spacecraft other than Cassini, you should provide
                the Allen deviation yourself. Default is sigma=2e-13
            :psitype (*str*):
                A string for determining what approximation to the
                geometrical 'psi' function is used. Several strings
                are allowed:

                |    'full'      No Approximation is applied.
                |    'MTR2'      Second Order Series from MTR86.
                |    'MTR3'      Third Order Series from MTR86.
                |    'MTR4'      Fourth Order Series from MTR86.
                |    'Fresnel'   Standard Fresnel approximation.

                The variable is neither case nor space sensitive.
                Default is set to 'full'.
            :verbose (*bool*):
                A Boolean for determining if various pieces of
                information are printed to the screen or not.
                Default is False.
        Attributes:
            :bfac (*bool*):
                Boolean for bfac (See keywords).
            :dathist (*dict*):
                History from DLP instance.
            :dx_km (*float*):
                Radial spacing for the data points (km).
            :f_sky_hz_vals (*np.ndarray*):
                Recieved frequency from the spacecraft (Hz).
            :finish (*int*):
                Final point that was reconstructed.
            :fwd (*bool*):
                Boolean for fwd (See keywords).
            :history (*dict*):
                History for the DiffractionCorrection class.
                This contains system info and user info, including
                what operating system was used, username, hostname,
                computer name, and the inputs provided.
            :lambda_sky_km_vals (*np.ndarray*):
                Wavelength of recieved signal from spacecraft (km).
            :mu_vals (*np.ndarray*):
                The sine of the ring opening angle (Unitless).
            :n_used (*int*):
                Number of points that were reconstructed.
            :norm (*bool*):
                Boolean for norm (See keywords).
            :norm_eq (*float*):
                Normalized equivalent width computed from window
                that was used during reconstruction. See the
                window_functions submodule for more information.
            :p_norm_fwd_vals (*np.ndarray*):
                Normalized power computer from the forward modelling
                of the reconstructed data. This will be a None type
                variable unless fwd=True is set. If the
                reconstruction went well, this should mimic the raw
                data, p_norm_vals.
            :p_norm_vals (*np.ndarray*):
                Normalized power from the diffracted signal. This is
                the square of the absolute value of the recieved
                complex transmittance.
            :phase_fwd_vals (*np.ndarray*):
                Phase computed from the forward model of the
                reconstructed data. This will be a None type
                variable unless fwd=True is set. If the
                reconstruction went well, this should mimic
                phase_rad_vals. This variable is in radians.
            :phase_rad_vals (*np.ndarray*):
                Phase from the diffracted signal (Radians).
            :phase_vals (*np.ndarray*):
                Reconstructed phase (Radians).
            :phi_rad_vals (*np.ndarray*):
                Ring azimuth angle of the ring intercept (Radians).
            :phi_rl_rad_vals (*np.ndarray*):
                Ring longitude angle. This will be a None type unless
                it was provided in the DLP class. Otherwise,
                this variable is in radians.
            :power_vals (*np.ndarray*):
                Normalized reconstructed power.
            :psitype (*str*):
                String for psitype (See keywords).
            :raw_tau_threshold_vals (*np.ndarray*):
                Threshold optical depth for the diffracted data.
                This will be a None type unless provided for in the
                DLP class.
            :res (*float*):
                Requested resolution (See arguments). In kilometers.
            :rho_corr_pole_km_vals (*np.ndarray*):
                Radial corrections from the Planet's pole. This will
                be a None type variable unless provided in the
                DLP class. Otherwise, this is in kilometers.
            :rho_corr_timing_km_vals (*np.ndarray*):
                Radial corrections from timing offsets. This will be
                a None type variable unless provided in the DLP
                class. Otherwise, this is in kilometers.
            :rho_dot_kms_vals (*np.ndarray*):
                Time derivative of the ring intercept point (km/s).
            :rho_km_vals (*np.ndarray*):
                Ring-intercept-point (RIP) in kilometers.
            :rng (*list*):
                Range that was used for reconstruction, taking into
                the range that was requested by the user. The actual
                range takes into account limits in the available data
                and limits in the required window sizes.
            :rngreq (*str* or *list*):
                Requested range (See keywords).
            :sigma (*float*):
                Requested Allen deviation (See keywords).
            :start (*int*):
                First point that was reconstructed.
            :t_oet_spm_vals (*np.ndarray*):
                Time the signal is measured on Earth. This is a
                None type unless provided for in the DLP class.
            :t_ret_spm_vals (*np.ndarray*):
                Time the signal passes through the diffracting
                medium. This is a None type unless provided for in
                the DLP class.
            :t_set_spm_vals (*np.ndarray*):
                Time the signal is emitted from the spacecraft. This
                is a None type unless provided in the DLP class.
            :tau_threshold_vals (*np.ndarray*):
                Threshold optical depth of the reconstructed data.
            :tau_vals (*np.ndarray*):
                Optical depth of the reconstructed data.
            :verbose (*bool*):
                Boolean for Verbose (See keywords).
            :w_km_vals (*np.ndarray*):
                Window width as a function of radius (km).
            :wtype (*str*):
                String for wtype (See keywords).
    """
    def __init__(self, DLP, res, rng="all", wtype="kbmd20", fwd=False,
                 norm=True, verbose=False, bfac=True, sigma=2.e-13,
                 psitype="fresnel4", write_file=False, res_factor=0.75,
                 eccentricity=0.0, periapse=0.0):

        fname = "diffrec.diffraction_correction.DiffractionCorrection"
        error_check.check_type(verbose, bool, "verbose", fname)

        if verbose:
            print("Processing Diffraction Correction:")
            print("\tRunning Error Check on Input Arguments...")
        else:
            pass

        error_check.check_type(fwd, bool, "fwd", fname)
        error_check.check_type(bfac, bool, "bfac", fname)
        error_check.check_type(norm, bool, "norm", fname)
        error_check.check_type(write_file, bool, "write_file", fname)

        res = error_check.check_type_and_convert(res, float, "res", fname)
        sigma = error_check.check_type_and_convert(sigma, float, "sigma", fname)
        res_factor = error_check.check_type_and_convert(res_factor, float,
                                                        "res_factor", fname)

        # Check that these variables are positive.
        error_check.check_positive(res, "res", fname)
        error_check.check_positive(sigma, "sigma", fname)
        error_check.check_positive(res_factor, "res_factor", fname)
        error_check.check_non_negative(eccentricity, "eccentricity", fname)
        
        # Check that the periapse is within [0, 2pi)
        error_check.check_two_pi(periapse, "periapse", fname)

        # Check that the requested window type is a legal input.
        if not isinstance(wtype, str):
            erm = ""
            for key in window_functions.func_dict:
                erm = "%s\t\t'%s'\n" % (erm, key)
            raise TypeError(
                "\n\tError Encountered: rss_ringoccs\n"
                "\t\tdiffrec.diffraction_correction.DiffractionCorrection\n\n"
                "\twtype must be a string.\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: str\n"
                "\tAllowed string are:\n%s" % (type(wtype).__name__, erm)
            )
        else:
            # Remove spaces and quotes from the wtype variable.
            wtype = wtype.replace(" ", "").replace("'", "").replace('"', "")

            # Set wtype string to lower-case.
            wtype = wtype.lower()
            if not (wtype in window_functions.func_dict):
                erm = ""
                for key in window_functions.func_dict:
                    erm = "%s\t\t'%s'\n" % (erm, key)
                raise ValueError(
                    "\n\tError Encountered: rss_ringoccs\n"
                    "\t\tdiffrec.diffraction_correction.DiffractionCorrection"
                    "\n\n\tIllegal string used for wtype.\n"
                    "\tYour string: '%s'\n"
                    "\tAllowed Strings:\n%s" % (wtype, erm)
                )
            else:
                pass

        # Check that range and psitype are legal inputs.
        rng = error_check.check_range_input(rng, fname)
        psitype = error_check.check_psitype(psitype, fname)

        if verbose:
            print("\tAssigning inputs as attributes...")

        # Assign variables as attributes.
        self.p_norm_fwd_vals = None
        self.T_hat_fwd_vals = None
        self.phase_fwd_vals = None
        self.eccentricity = eccentricity
        self.periapse = periapse
        self.input_res = res
        self.verbose = verbose
        self.psitype = psitype
        self.rngreq = rng
        self.wtype = wtype
        self.sigma = sigma
        self.norm = norm
        self.bfac = bfac
        self.res = res*res_factor
        self.fwd = fwd

        # Retrieve variables from the DLP class, setting as attributes.
        if verbose:
            print("\tRetrieving variables from DLP instance...")

        try:
            erm = "rho_km_vals"
            self.rho_km_vals = np.array(DLP.rho_km_vals)
            erm = "p_norm_vals"
            self.p_norm_vals = np.array(DLP.p_norm_vals)
            erm = "phase_rad_vals"
            self.phase_rad_vals = np.array(DLP.phase_rad_vals)
            erm = "B_rad_vals"
            self.B_rad_vals = np.array(DLP.B_rad_vals)
            erm = "D_km_vals"
            self.D_km_vals = np.array(DLP.D_km_vals)
            erm = "phi_rad_vals"
            self.phi_rad_vals = np.array(DLP.phi_rad_vals)
            erm = "f_sky_hz_vals"
            self.f_sky_hz_vals = np.array(DLP.f_sky_hz_vals)
            erm = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(DLP.rho_dot_kms_vals)
            erm = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(DLP.t_oet_spm_vals)
            erm = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(DLP.t_ret_spm_vals)
            erm = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(DLP.t_set_spm_vals)
            erm = "rho_corr_pole_km_vals"
            self.rho_corr_pole_km_vals = np.array(DLP.rho_corr_pole_km_vals)
            erm = "rho_corr_timing_km_vals"
            self.rho_corr_timing_km_vals = np.array(DLP.rho_corr_timing_km_vals)
            erm = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = np.array(DLP.phi_rl_rad_vals)
            erm = "raw_tau_threshold_vals"
            self.raw_tau_threshold_vals = np.array(DLP.raw_tau_threshold_vals)
        except (TypeError, ValueError, NameError, AttributeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\t\trss_ringoccs.diffrec.DiffractionCorrection\n\n"
                "\t%s could not be converted into a numpy array.\n"
                "\tCheck your DLP class for errors." % erm
            )
        
        self.dathist = DLP.history

        # Run various error checks on all variables.
        error_check.check_is_real(self.D_km_vals, "D_km_vals", fname)
        error_check.check_is_real(self.B_rad_vals, "B_rad_vals", fname)
        error_check.check_is_real(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_is_real(self.p_norm_vals, "p_norm_vals", fname)
        error_check.check_is_real(self.phi_rad_vals, "phi_rad_vals", fname)
        error_check.check_is_real(self.f_sky_hz_vals, "f_sky_hz_vals", fname)
        error_check.check_is_real(self.phase_rad_vals, "phase_rad_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals, 
                                  "rho_dot_kms_vals", fname)

        error_check.check_positive(self.D_km_vals, "D_km_vals", fname)
        error_check.check_positive(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_positive(self.f_sky_hz_vals, "f_sky_hz_vals", fname)
        error_check.check_non_negative(self.p_norm_vals, "p_norm_vals", fname)

        error_check.check_two_pi(self.B_rad_vals, "B_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rad_vals, "phi_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phase_rad_vals, "phase_rad_vals",
                                 fname, deg=False)

        if (np.size(self.rho_km_vals) < 2):
            raise IndexError(
                """
                    \r\tError Encountered:
                    \r\t\trss_ringoccs.diffrec.DiffractionCorrection\n
                    \r\trho_km_vals has less than 2 points.
                    \r\tIt is impossible to do reconstruction.
                """
            )
        else:
            self.rho_km_vals = self.rho_km_vals.astype(float)

        error_check.check_lengths(self.rho_dot_kms_vals, self.rho_km_vals, 
                                  "rho_dot_kms_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phase_rad_vals, self.rho_km_vals,
                                  "phase_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.f_sky_hz_vals, self.rho_km_vals,
                                  "f_sky_hz_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.phi_rad_vals, self.rho_km_vals,
                                  "phi_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.p_norm_vals, self.rho_km_vals,
                                  "p_norm_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.B_rad_vals, self.rho_km_vals,
                                  "B_rad_vals", "rho_km_vals", fname)
        error_check.check_lengths(self.D_km_vals, self.rho_km_vals,
                                  "D_km_vals", "rho_km_vals", fname)
        
        self.rho_dot_kms_vals = self.rho_dot_kms_vals.astype(float)
        self.phase_rad_vals = -self.phase_rad_vals.astype(float)
        self.f_sky_hz_vals = self.f_sky_hz_vals.astype(float)
        self.phi_rad_vals = self.phi_rad_vals.astype(float)
        self.p_norm_vals = self.p_norm_vals.astype(float)
        self.B_rad_vals = self.B_rad_vals.astype(float)
        self.D_km_vals = self.D_km_vals.astype(float)

        # Compute sampling distance (km)
        self.dx_km = self.rho_km_vals[1] - self.rho_km_vals[0]

        # Check that the data is well sampled for the requested resolution.
        if (self.dx_km == 0.0):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\trho_km_vals[1]-rho_km_vals[0]=0.0
                    \r\tThe sample spacing is zero.
                """ % (fname)
            )
        elif self.res < 1.999999*self.dx_km:
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tRequested resolution is less than twice the sample
                    \r\tspacing. This will produce an inaccurate results.\n
                    \r\tRequested Resolution (km): %f
                    \r\tSample Spacing (km): %f\n
                    \r\tTO CORRECT THIS:
                    \r\t\tChoose a resolution GREATER than %f km\n
                    \r\tPLEASE NOTE:\n
                    \r\t\tTo be consistent with PDS results, a scale factor
                    \r\t\tof 0.75 is applied to your requested resolution.
                    \r\t\tto ignore this, set 'res_factor=1.0' when calling
                    \r\t\tthe DiffractionCorrection class.\n
                    \r\t\tres_factor is currently set to: %f
                """ % (fname, self.res, self.dx_km,
                       2.0*self.dx_km/res_factor, res_factor)
            )
        else:
            pass

        if verbose:
            print("\tCheck Variables for Errors...")

        # Check that rho_km_vals is increasing and the rev isn't a chord occ.
        drho = [np.min(self.rho_dot_kms_vals), np.max(self.rho_dot_kms_vals)]

        if (drho[0] < 0) and (drho[1] > 0):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tdrho/dt has positive and negative values.
                    \r\tYour input file is probably a chord occultation.
                    \r\tDiffraction Correction can only be performed for
                    \r\tone event at a time. That is, ingress or egress.\n
                    \r\tTO CORRECT THIS:
                    \r\t\tSplit the input into two parts: An egress
                    \r\t\tportion and an ingress portion, and then run
                    \r\t\tdiffraction correction on the individual pieces.
                """ % (fname)
            )
        elif ((drho[0] == 0.0) or (drho[1] == 0.0)):
            raise ValueError(
                """
                    \r\tError Encountered:
                    \r\t\t%s\n
                    \r\tdrho/dt has elements with value zero.
                    \r\tYour input file is probably a chord occultation.
                    \r\tDiffraction Correction can only be performed for
                    \r\tone event at a time. That is, either an ingress
                    \r\tor an egress event.\n
                    \r\tTO CORRECT THIS:
                    \r\t\tSplit the input into two parts: An egress
                    \r\t\tportion and an ingress portion, and then run
                    \r\t\tdiffraction correction on the individual pieces.
                    \r\t\tIgnore the region where drho/dt is close to zero.
                """ % (fname)
            )
        elif (self.dx_km > 0) and (drho[1] < 0):
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals)
        elif (self.dx_km < 0) and (drho[0] > 0):
            raise ValueError(
                """
                    \r\tError Encountered:
                    \r\t\t%s\n
                    \r\trho_km_vals is decreasing yet rho_dot_kms_vals
                    \r\tis positiive. Check DLP class for errors.
                """ % (fname)
            )
        elif (self.dx_km < 0):
            self.rho_km_vals = self.rho_km_vals[::-1]
            self.phase_rad_vals = self.phase_rad_vals[::-1]
            self.p_norm_vals = self.p_norm_vals[::-1]
            self.phi_rad_vals = self.phi_rad_vals[::-1]
            self.B_rad_vals = self.B_rad_vals[::-1]
            self.f_sky_hz_vals = self.f_sky_hz_vals[::-1]
            self.D_km_vals = self.D_km_vals[::-1]
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])
            self.dx_km *= -1.0
        else:
            del drho

        if verbose:
            print("\tComputing Necessary Variables...")

        # Compute various variables.
        self.lambda_sky_km_vals = SPEED_OF_LIGHT_KM / self.f_sky_hz_vals
        self.mu_vals = np.sin(np.abs(self.B_rad_vals))
        self.T_hat_vals = np.exp(1j*self.phase_rad_vals)
        self.T_hat_vals *= np.sqrt(self.p_norm_vals)

        # Compute geometric qunatities and the Fresnel Scale.
        self.F_km_vals = special_functions.fresnel_scale(
            self.lambda_sky_km_vals, self.D_km_vals,
            self.phi_rad_vals, self.B_rad_vals, deg=False
        )

        # Compute the Normalized Equaivalent Width (See MTR86 Equation 20)
        self.norm_eq = window_functions.func_dict[wtype]["normeq"]

        # Compute the window width. (See MTR86 Equations 19, 32, and 33).
        self.w_km_vals, Prange = window_functions.window_width(
            self.res, self.norm_eq, self.f_sky_hz_vals, self.F_km_vals,
            self.rho_dot_kms_vals, self.sigma, bfac=self.bfac, Return_P=True
        )

        # From the requested range, extract array of the form [a, b]
        if (isinstance(rng, str)):
            self.rng = np.array(region_dict[rng])
        else:
            self.rng = np.array([np.min(rng), np.max(rng)])

        # Compute the smallest and largest allowed radii for reconstruction.
        rho = self.rho_km_vals[Prange]
        w = self.w_km_vals[Prange]
        rho_min = self.rho_km_vals[Prange]-self.w_km_vals[Prange]/2.0
        rho_max = self.rho_km_vals[Prange]+self.w_km_vals[Prange]/2.0

        wrange = Prange[np.where((rho_min >= np.min(rho)) &
                                 (rho_max <= np.max(rho)))]
        self.wrange = wrange

        # Check that there is enough data for reconstruction.
        if (np.size(wrange) == 0):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tThe window width is too large to reconstruct any
                    \r\tpoints. Please choose a coarser resolution or
                    \r\tinspect your input data.\n
                    \r\t\tMinimum Available Radius:         %f
                    \r\t\tMaximum Available Radius:         %f
                    \r\t\tMinimum Required Window Width:    %f
                    \r\t\tMaximum Required Window Width:    %f
                """ % (fname, np.min(rho), np.max(rho), np.min(w), np.max(w))
            )
        elif (np.max(rho) < np.min(self.rng)):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tMinimum requested range is greater
                    \r\tthan the maximum available data point.
                    \r\tSelect a smaller range for reconstruction.\n
                    \r\tYour Requested Minimum (km):    %f
                    \r\tYour Requested Maximum (km):    %f
                    \r\tMaximum Available Data (km):    %f
                """ % (fname, np.min(self.rng), np.max(self.rng), np.max(rho))
            )
        elif (np.min(rho) > np.max(self.rng)):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tMaximum requested range is less
                    \r\tthan the minimum available data point.\n
                    \r\tYour Requested Minimum (km): %f
                    \r\tYour Requested Maximum (km): %f
                    \r\tMinimum Available Data (km): %f\n
                    \r\tTO CORRECT THIS:
                    \r\t\tSelect a larger range for reconstruction
                """ % (fname, np.min(self.rng), np.max(self.rng), np.min(rho))
            )
        else:
            pass

        rho_min = np.min(rho[wrange])
        rho_max = np.max(rho[wrange])

        wrange = wrange[np.where((rho[wrange] >= np.min(self.rng)) &
                                 (rho[wrange] <= np.max(self.rng)))]

        if (np.size(wrange) <= 1):
            raise IndexError(
                "\n\tError Encountered:\n"
                "\t\trss_ringoccs.diffrec.DiffractionCorrection\n\n"
                "\tRequested range does not include any of the\n"
                "\tavailable points for processing. Please choose\n"
                "\tA different range for processing.\n"
                "\t\tMinimum Possible Radius: %f\n"
                "\t\tMaximum Possible Radius: %f\n"
                "\t\tRequested Range Minimum: %f\n"
                "\t\tRequested Range Maximum: %f"
                % (rho_min, rho_max, np.min(self.rng), np.max(self.rng))
            )
        else:
            pass

        self.start = wrange[0]
        self.finish = wrange[-1]
        self.n_used = 1 + (self.finish - self.start)

        # Create input variable and keyword dictionaries for history.
        input_vars = {
            'dlp_inst': DLP.history,
            'res':      res
        }

        input_kwds = {
            'rng':          rng,
            'wtype':        wtype,
            'fwd':          fwd,
            'norm':         norm,
            'bfac':         bfac,
            'sigma':        sigma,
            'psitype':      psitype,
            'res_factor':   res_factor,
            'periapse':     periapse,
            'eccentricity': eccentricity
        }

        # Delete unnecessary variables for clarity.
        del rho, rho_min, rho_max, rng, norm, fwd, bfac, psitype, verbose

        if self.verbose:
            print("\tRunning Fresnel Inversion...")

        self.T_vals = self.__ftrans(fwd=False)

        # Compute power and phase.
        if self.verbose:
            print("\tComputing Power and Phase...")

        self.power_vals = np.square(np.abs(self.T_vals))

        # Return phase to original sign.
        self.phase_vals = -np.arctan2(np.imag(self.T_vals),
                                      np.real(self.T_vals))
        self.phase_rad_vals *= -1

        if self.verbose:
            print("\tInversion Complete.")

        if self.fwd:
            if self.verbose:
                print("\tComputing Forward Transform...")

            self.T_hat_fwd_vals = self.__ftrans(fwd=True)
            self.p_norm_fwd_vals = np.square(np.abs(self.T_hat_fwd_vals))
            self.phase_fwd_vals = -np.arctan2(np.imag(self.T_hat_fwd_vals),
                                              np.real(self.T_hat_fwd_vals))
            if self.verbose:
                print("\tForward Transform Complete.")

        # Compute regions of non-zero power.
        crange = (self.power_vals > 0.0).nonzero()

        # Compute the normalized optical depth.
        self.tau_vals = np.zeros(np.size(self.power_vals))
        self.tau_vals[crange] = np.log(self.power_vals[crange])
        self.tau_vals[crange] *= -self.mu_vals[crange]

        self.tau_threshold_vals = (self.raw_tau_threshold_vals -
                                   self.mu_vals*np.log(self.dx_km/self.res))

        self.__trim_attributes(self.fwd)

        self.history = write_history_dict(input_vars, input_kwds, __file__)

        # Set rev_info attribute from DLP instance.
        self.rev_info = DLP.rev_info
        if write_file:
            self.outfiles = write_output_files(self)

        if self.verbose:
            print("\tDiffraction Correction Complete.")

    def __trim_attributes(self, fwd):
        """
            Purpose:
                Trim the attributes in the DiffractionCorrection
                class so that only reconstructed points will be
                returned to the user. All other unused points are
                discarded.
            Keywords:
                :fwd (*bool*):
                    Boolean for the forward calculation.
                    If set to True, the forward variables
                    will also be trimmed.
        """
        # Get rid of uncomputed values and keep only what was processed.
        start = self.start
        n_used = self.n_used
        crange = np.arange(n_used)+start

        # Ring radius, azimuth angle, diffracted power, and phase.
        self.rho_km_vals = self.rho_km_vals[crange]
        self.p_norm_vals = self.p_norm_vals[crange]
        self.phi_rad_vals = self.phi_rad_vals[crange]
        self.phase_rad_vals = self.phase_rad_vals[crange]

        # Ring opening angle, normalized power, phase, and transmittance.
        self.B_rad_vals = self.B_rad_vals[crange]
        self.power_vals = self.power_vals[crange]
        self.phase_vals = self.phase_vals[crange]
        self.T_hat_vals = self.T_hat_vals[crange]
        self.T_vals = self.T_vals[crange]

        # Fresnel scale, window width, and RIP-Spacecraft distance.
        self.F_km_vals = self.F_km_vals[crange]
        self.w_km_vals = self.w_km_vals[crange]
        self.D_km_vals = self.D_km_vals[crange]

        # Ring radius corrections.
        self.rho_corr_pole_km_vals = self.rho_corr_pole_km_vals[crange]
        self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals[crange]

        # Various time attributes.
        self.t_oet_spm_vals = self.t_oet_spm_vals[crange]
        self.t_ret_spm_vals = self.t_ret_spm_vals[crange]
        self.t_set_spm_vals = self.t_set_spm_vals[crange]

        # All other attributes.
        self.mu_vals = self.mu_vals[crange]
        self.tau_vals = self.tau_vals[crange]
        self.f_sky_hz_vals = self.f_sky_hz_vals[crange]
        self.phi_rl_rad_vals = self.phi_rl_rad_vals[crange]
        self.rho_dot_kms_vals = self.rho_dot_kms_vals[crange]
        self.lambda_sky_km_vals = self.lambda_sky_km_vals[crange]
        self.raw_tau_threshold_vals = self.raw_tau_threshold_vals[crange]
        self.tau_threshold_vals = self.tau_threshold_vals[crange]

        # If the forward model was run, trim those attributes as well.
        if fwd:
            # Forward power
            self.p_norm_fwd_vals = self.p_norm_fwd_vals[crange]

            # Forward Transmittance and phase.
            self.T_hat_fwd_vals = self.T_hat_fwd_vals[crange]
            self.phase_fwd_vals = self.phase_fwd_vals[crange]

    def __ftrans(self, fwd):
        """
            Purpose:
                Compute the Fresnel Inversion.
            Arguments:
                :self:
                    Instance of DiffractionCorrection class.
            Keywords:
                :fwd (*bool*):
                    Boolean for whether or not the forward
                    calculation is being performed.
            Outputs:
                :T_out (*np.ndarray*):
                    Complex transmittance.
        """
        # Compute product of wavenumber and RIP distance.
        kD_vals = TWO_PI * self.D_km_vals / self.lambda_sky_km_vals

        # Define functions.
        fw = window_functions.func_dict[self.wtype]["func"]
        mes = "\t\tPt: %d  Tot: %d  Width: %d  Psi Iters: %d\t"
        __norm = window_functions.normalize

        # If forward transform, adjust starting point by half a window.
        if fwd:
            w_max = np.max(self.w_km_vals[self.start:self.start + self.n_used])
            nw_fwd = int(np.ceil(w_max / (2.0 * self.dx_km)))
            start = int(self.start + nw_fwd)
            n_used = int(self.n_used - 2 * nw_fwd)
            T_in = self.T_vals
        else:
            start = self.start
            n_used = self.n_used
            T_in = self.T_hat_vals
        
        # Create empty array for reconstruction / forward transform.
        T_out = T_in * 0.0

        # Compute first window width and window function.
        w_init = self.w_km_vals[start]
        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
        crange = np.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))
        r0 = self.rho_km_vals[crange]
        r = self.rho_km_vals[start]
        x = r-r0
        w_func = fw(x, w_init)

        # Assign integer to window functions for passing args into C code.
        if (self.wtype == "rect"):
            wnum = 0
        elif (self.wtype == "coss"):
            wnum = 1
        elif (self.wtype == "kb20"):
            wnum = 2
        elif (self.wtype == "kb25"):
            wnum = 3
        elif (self.wtype == "kb35"):
            wnum = 4
        elif (self.wtype == "kbmd20"):
            wnum = 5
        elif (self.wtype == "kbmd25"):
            wnum = 6
        else:
            wnum = 7

        if (self.norm == False):
            use_norm = 0
        else:
            use_norm = 1

        if (fwd == False):
            use_fwd = 0
        else:
            use_fwd = 1

        if (self.psitype == "ellipse"):
            peri = self.periapse
            ecc = self.eccentricity
            for i in np.arange(n_used):
                # Current point being computed.
                center = start+i

                # Current window width, Fresnel scale, and ring radius.
                w = self.w_km_vals[center]
                F = self.F_km_vals[center]
                r = self.rho_km_vals[center]

                if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                    # Compute first window width and window function.
                    w_init = self.w_km_vals[center]

                    nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
                    crange = np.arange(int(center-(nw-1)/2),
                                       int(1+center+(nw-1)/2))

                    # Ajdust ring radius by dx_km.
                    r0 = self.rho_km_vals[crange]

                    # Compute psi for with stationary phase value
                    x = r-r0
                    w_func = fw(x, w_init)
                else:
                    crange += 1
                    r0 = self.rho_km_vals[crange]
                
                d = self.D_km_vals[center]
                b = self.B_rad_vals[center]
                kD = kD_vals[center]
                phi = self.phi_rad_vals[crange]
                phi0 = self.phi_rad_vals[crange]

                # Compute Newton-Raphson perturbation
                psi_d1 = special_functions.dpsi_ellipse(kD, r, r0, phi, phi0,
                                                        b, d, ecc, peri)
                loop = 0

                while (np.max(np.abs(psi_d1)) > 1.0e-4):
                    psi_d1 = special_functions.dpsi_ellipse(kD, r, r0,
                                                            phi, phi0, b,
                                                            d, ecc, peri)
                    psi_d2 = special_functions.d2psi(kD, r, r0, phi, phi0, b, d)
                    
                    # Newton-Raphson
                    phi += -(psi_d1 / psi_d2)

                    # Add one to loop variable for each iteration
                    loop += 1
                    if (loop > 4):
                        break

                # Compute Eta variable (MTR86 Equation 4c).
                psi_vals = special_functions.psi(kD, r, r0, phi, phi0, b, d)

                # Compute kernel function for Fresnel inverse
                if fwd:
                    ker = w_func*np.exp(1j*psi_vals)
                else:
                    ker = w_func*np.exp(-1j*psi_vals)

                # Compute 'approximate' Fresnel Inversion for current point
                T = T_in[crange]
                T_out[center] = np.sum(ker*T)*self.dx_km*(0.5+0.5j)/F

                # If normalization has been set, normalize the reconstruction
                if self.norm:
                    T_out[center] *= window_functions.normalize(self.dx_km, ker, F)
                if self.verbose:
                    print(mes % (i, n_used-1, nw, loop), end="\r")
        elif (self.psitype == "full"):
            try:
                if self.verbose:
                    print("\t\tTrying C Code for Diffraction Correction...")

                return _diffraction_functions.fresnel_transform_newton(
                    T_in, self.rho_km_vals, self.F_km_vals,
                    self.phi_rad_vals, kD_vals, self.B_rad_vals,
                    self.D_km_vals, self.w_km_vals, start, n_used,
                    wnum, use_norm, use_fwd
                )
            except:
                for i in np.arange(n_used):
                    # Current point being computed.
                    center = start+i

                    # Current window width, Fresnel scale, and ring radius.
                    w = self.w_km_vals[center]
                    F = self.F_km_vals[center]
                    r = self.rho_km_vals[center]

                    if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                        # Compute first window width and window function.
                        w_init = self.w_km_vals[center]

                        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
                        crange = np.arange(int(center-(nw-1)/2),
                                        int(1+center+(nw-1)/2))

                        # Ajdust ring radius by dx_km.
                        r0 = self.rho_km_vals[crange]

                        # Compute psi for with stationary phase value
                        x = r-r0
                        w_func = fw(x, w_init)
                    else:
                        crange += 1
                        r0 = self.rho_km_vals[crange]
                    
                    d = self.D_km_vals[center]
                    b = self.B_rad_vals[center]
                    kD = kD_vals[center]
                    phi = self.phi_rad_vals[center]
                    phi0 = self.phi_rad_vals[center]

                    # Compute Newton-Raphson perturbation
                    psi_d1 = special_functions.dpsi(kD, r, r0, phi, phi0, b, d)
                    loop = 0

                    while (np.max(np.abs(psi_d1)) > 1.0e-4):
                        psi_d1 = special_functions.dpsi(kD, r, r0, phi, phi0, b, d)
                        psi_d2 = special_functions.d2psi(kD, r, r0, phi, phi0, b, d)
                        
                        # Newton-Raphson
                        phi += -(psi_d1 / psi_d2)

                        # Add one to loop variable for each iteration
                        loop += 1
                        if (loop > 4):
                            break

                    # Compute Eta variable (MTR86 Equation 4c).
                    psi_vals = special_functions.psi(kD, r, r0, phi, phi0, b, d)

                    # Compute kernel function for Fresnel inverse
                    if fwd:
                        ker = w_func*np.exp(1j*psi_vals)
                    else:
                        ker = w_func*np.exp(-1j*psi_vals)

                    # Compute 'approximate' Fresnel Inversion for current point
                    T = T_in[crange]
                    T_out[center] = np.sum(ker*T)*self.dx_km*(0.5+0.5j)/F

                    # If normalization has been set, normalize the reconstruction
                    if self.norm:
                        T_out[center] *= window_functions.normalize(self.dx_km, ker, F)
                    if self.verbose:
                        print(mes % (i, n_used-1, nw, loop), end="\r")
        else:
            if (self.psitype == "fresnel"):
                try:
                    if self.verbose:
                        print("\t\tTrying C Code for Diffraction Correction...")

                    return _diffraction_functions.fresnel_transform_quadratic(
                        T_in, self.dx_km, self.F_km_vals, self.w_km_vals,
                        start, n_used, wnum, use_norm, use_fwd
                    )
                except KeyboardInterrupt:
                    sys.exit("KeyboardInterrupt")
                except:
                    if self.verbose:
                        print("\t\tCould not import C code. Using Python Code.")

                    crange -= 1
                    F2 = np.square(self.F_km_vals)
                    x = r-r0
                    x2 = HALF_PI * np.square(x)
                    loop = 0
                    for i in np.arange(n_used):
                        # Current point being computed.
                        center = start+i

                        # Window width and Frensel scale for current point.
                        w = self.w_km_vals[center]
                        F = self.F_km_vals[center]

                        if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                            # Compute first window width and window function.
                            w_init = self.w_km_vals[center]
                            nw = int(2 * np.floor(w_init / (2.0*self.dx_km))+1)        
                            crange = np.arange(int(center-(nw-1)/2),
                                            int(1+center+(nw-1)/2))

                            # Ajdust ring radius by dx_km.
                            r0 = self.rho_km_vals[crange]
                            r = self.rho_km_vals[center]

                            # Compute psi for with stationary phase value
                            x = r-r0
                            x2 = HALF_PI * np.square(x)

                            w_func = fw(x, w_init)
                        else:
                            crange += 1
                        
                        psi_vals = x2 / F2[center]

                        # Compute kernel function for Fresnel inverse
                        if fwd:
                            ker = w_func*np.exp(1j*psi_vals)
                        else:
                            ker = w_func*np.exp(-1j*psi_vals)

                        # Range of diffracted data that falls inside the window
                        T = T_in[crange]

                        # Compute 'approximate' Fresnel Inversion.
                        T_out[center] = np.sum(ker*T)*self.dx_km*(0.5+0.5j)/F

                        # If norm has been set, normalize the reconstruction
                        if self.norm:
                            T_out[center] *= __norm(self.dx_km, ker, F)
                        if self.verbose:
                            print(mes % (i, n_used, nw, loop), end="\r")
                    if self.verbose:
                        print(mes % (i, n_used, nw, loop))
                    return T_out
            elif (self.psitype == "fresnel3"):
                try:
                    if self.verbose:
                        print("\t\tTrying C Code for Diffraction Correction...")
                    return _diffraction_functions.fresnel_transform_cubic(
                        T_in, self.dx_km, self.F_km_vals, self.phi_rad_vals,
                        kD_vals, self.B_rad_vals, self.D_km_vals,
                        self.w_km_vals, start, n_used, wnum, use_norm, use_fwd
                    )
                except (TypeError, ValueError, NameError):
                    if self.verbose:
                        print("\t\tCould not import C code. Using Python Code.")
            elif (self.psitype == "fresnel4"):
                try:
                    if self.verbose:
                        print("\t\tTrying C Code for Diffraction Correction...")

                    return _diffraction_functions.fresnel_transform_quartic(
                        T_in, self.dx_km, self.F_km_vals, self.phi_rad_vals,
                        kD_vals, self.B_rad_vals, self.D_km_vals,
                        self.w_km_vals, start, n_used, wnum, use_norm, use_fwd
                    )
                except (TypeError, ValueError, NameError):   
                    if self.verbose:
                        print("\t\tCould not import C code. Using Python Code.")
            elif (self.psitype == "fresnel6"):
                try:
                    if self.verbose:
                        print("\t\tTrying C Code for Diffraction Correction...")

                    return _diffraction_functions.fresnel_transform_sextic(
                        T_in, self.dx_km, self.F_km_vals, self.phi_rad_vals,
                        kD_vals, self.B_rad_vals, self.D_km_vals,
                        self.w_km_vals, start, n_used, wnum, use_norm, use_fwd
                    )
                except (TypeError, ValueError, NameError):
                    if self.verbose:
                        print("\t\tCould not import C code. Using Python Code.")
            else:
                try:
                    if self.verbose:
                        print("\t\tTrying C Code for Diffraction Correction...")

                    return _diffraction_functions.fresnel_transform_octic(
                        T_in, self.dx_km, self.F_km_vals, self.phi_rad_vals,
                        kD_vals, self.B_rad_vals, self.D_km_vals,
                        self.w_km_vals, start, n_used, wnum, use_norm, use_fwd
                    )
                except (TypeError, ValueError, NameError):
                    if self.verbose:
                        print("\t\tCould not import C code. Using Python Code.")

            crange -= 1
            cosb = np.cos(self.B_rad_vals)
            cosp = np.cos(self.phi_rad_vals)
            sinp = np.sin(self.phi_rad_vals)
            A_2 = 0.5*np.square(cosb*sinp)/(1.0-np.square(cosb*sinp))

            # Legendre Polynomials.
            P_1 = cosb*cosp;
            P12 = P_1*P_1;
            P_2 = 1.5*P12-0.5;
            P_3 = (2.5*P12-1.5)*P_1;
            P_4 = (4.375*P12-7.5)*P12+0.376;
            P_5 = ((7.875*P12-8.75)*P12+1.875)*P_1;
            P_6 = (((14.4375*P12-19.6875)*P12+6.5625)*P12)-0.3125;
            P_7 = (((P12*26.8125-43.3125)*P12+19.6875)*P12-2.1875)*P_1;

            # Second set of polynomials.
            b_0 = 0.5 - 0.5*P12
            b_1 = (0.333333333333 - 0.333333333333*P_2)*P_1
            b_2 = (P_2 - P_1*P_3)*0.25
            b_3 = (P_3 - P_1*P_4)*0.20
            b_4 = (P_4 - P_1*P_5)*0.16666666666666666
            b_5 = (P_5 - P_1*P_6)*0.14285714285714285
            b_6 = (P_6 - P_1*P_7)*0.125

            # Products of Legendre Polynomials used in Expansion
            C_0 = P12
            C_1 = 2.0*P_1*P_2
            C_2 = 2.0*P_1*P_3+np.square(P_2)
            C_3 = 2.0*P_1*P_4+2.0*P_2*P_3
            C_4 = 2.0*P_2*P_4+np.square(P_3)
            C_5 = 2.0*P_3*P_4
            C_6 = np.square(P_4)

            x = (r-r0)
            x2 = np.square(x)
            x3 = x2*x
            x4 = np.square(x2)
            x5 = x4*x
            x6 = np.square(x3)

            # D_km_vals and various powers.
            d = self.D_km_vals
            d2 = np.square(d)
            d3 = d2*d
            d4 = np.square(d2)
            d5 = d4*d
            d6 = np.square(d3)

            loop = 0
            if (self.psitype == "fresnel3"):
                for i in np.arange(n_used):
                    # Current point being computed.
                    center = start+i

                    # Window width and Frensel scale for current point.
                    w = self.w_km_vals[center]
                    F = self.F_km_vals[center]

                    if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                        # Compute first window width and window function.
                        w_init = self.w_km_vals[center]
                        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
                        crange = np.arange(int(center-(nw-1)/2),
                                           int(1+center+(nw-1)/2))

                        # Ajdust ring radius by dx_km.
                        r0 = self.rho_km_vals[crange]
                        r = self.rho_km_vals[center]

                        # Compute psi for with stationary phase value
                        x = r-r0
                        x2 = np.square(x)

                        w_func = fw(x, w_init)
                    else:
                        crange += 1
                    
                    z = x/d[center]
                    z2 = x2/d2[center]

                    psi_vals = z2*(b_0[center]-A_2[center]*C_0[center]+
                                  (b_1[center]-A_2[center]*C_1[center])*z)
                    psi_vals *= kD_vals[center]

                    # Compute kernel function for Fresnel inverse
                    if fwd:
                        ker = w_func*np.exp(1j*psi_vals)
                    else:
                        ker = w_func*np.exp(-1j*psi_vals)

                    # Range of diffracted data that falls inside the window
                    T = T_in[crange]

                    # Compute 'approximate' Fresnel Inversion for current point
                    T_out[center] = np.sum(ker*T)*self.dx_km*(1.0+1.0j)/(2.0*F)

                    # If normalization has been set, normalize correction.
                    if self.norm:
                        T_out[center] *= __norm(self.dx_km, ker, F)
                    if self.verbose:
                        print(mes % (i, n_used, nw, loop), end="\r")
            elif (self.psitype == "fresnel4"):
                for i in np.arange(n_used):
                    # Current point being computed.
                    center = start+i

                    # Window width and Frensel scale for current point.
                    w = self.w_km_vals[center]
                    F = self.F_km_vals[center]

                    if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                        # Compute first window width and window function.
                        w_init = self.w_km_vals[center]
                        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
                        crange = np.arange(int(center-(nw-1)/2),
                                        int(1+center+(nw-1)/2))

                        # Ajdust ring radius by dx_km.
                        r0 = self.rho_km_vals[crange]
                        r = self.rho_km_vals[center]

                        # Compute psi for with stationary phase value
                        x = r-r0
                        x2 = np.square(x)

                        w_func = fw(x, w_init)
                    else:
                        crange += 1

                    z = x/d[center]
                    z2 = x2/d2[center]

                    psi_vals = b_0[center]-A_2[center]*C_0[center]
                    psi_vals += z*(b_1[center]-A_2[center]*C_1[center])
                    psi_vals += z2*(b_2[center]-A_2[center]*C_2[center])
                    psi_vals *= kD_vals[center]*z2

                    # Compute kernel function for Fresnel inverse
                    if fwd:
                        ker = w_func*np.exp(1j*psi_vals)
                    else:
                        ker = w_func*np.exp(-1j*psi_vals)

                    # Range of diffracted data that falls inside the window
                    T = T_in[crange]

                    # Compute 'approximate' Fresnel Inversion for current point
                    T_out[center] = np.sum(ker*T)*self.dx_km*(1.0+1.0j)/(2.0*F)

                    # If normalization has been set, normalize the reconstruction
                    if self.norm:
                        T_out[center] *= __norm(self.dx_km, ker, F)
                    if self.verbose:
                        print(mes % (i, n_used, nw, loop), end="\r")
            elif (self.psitype == "fresnel6"):
                for i in np.arange(n_used):
                    # Current point being computed.
                    center = start+i

                    # Window width and Frensel scale for current point.
                    w = self.w_km_vals[center]
                    F = self.F_km_vals[center]

                    if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                        # Compute first window width and window function.
                        w_init = self.w_km_vals[center]
                        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)        
                        crange = np.arange(int(center-(nw-1)/2),
                                        int(1+center+(nw-1)/2))

                        # Ajdust ring radius by dx_km.
                        r0 = self.rho_km_vals[crange]
                        r = self.rho_km_vals[center]

                        # Compute psi for with stationary phase value
                        x = r-r0
                        x2 = np.square(x)
                        x3 = x2*x
                        x4 = np.square(x2)

                        w_func = fw(x, w_init)
                    else:
                        crange += 1

                    z = x/d[center]
                    z2 = x2/d2[center]
                    z3 = x3/d3[center]
                    z4 = x4/d4[center]

                    psi_vals = b_0[center]-A_2[crange]*C_0[center]
                    psi_vals += z*(b_1[center]-A_2[crange]*C_1[center])
                    psi_vals += z2*(b_2[center]-A_2[crange]*C_2[center])
                    psi_vals += z3*(b_3[center]-A_2[crange]*C_3[center])
                    psi_vals += z4*(b_4[center]-A_2[crange]*C_4[center])

                    psi_vals *= z2*kD_vals[center]

                    # Compute kernel function for Fresnel inverse
                    if fwd:
                        ker = w_func*np.exp(1j*psi_vals)
                    else:
                        ker = w_func*np.exp(-1j*psi_vals)

                    # Range of diffracted data that falls inside the window
                    T = T_in[crange]

                    # Compute 'approximate' Fresnel Inversion for current point
                    T_out[center] = np.sum(ker*T)*self.dx_km*(1.0+1.0j)/(2.0*F)

                    # If normalization has been set, normalize the reconstruction
                    if self.norm:
                        T_out[center] *= __norm(self.dx_km, ker, F)
                    if self.verbose:
                        print(mes % (i, n_used, nw, loop), end="\r")
            else:
                for i in np.arange(n_used):
                    # Current point being computed.
                    center = start+i

                    # Window width and Frensel scale for current point.
                    w = self.w_km_vals[center]
                    F = self.F_km_vals[center]

                    if (np.abs(w_init - w) >= 2.0 * self.dx_km):

                        # Compute first window width and window function.
                        w_init = self.w_km_vals[center]
                        nw = int(2 * np.floor(w_init / (2.0 * self.dx_km)) + 1)
                        crange = np.arange(int(center-(nw-1)/2),
                                           int(1+center+(nw-1)/2))

                        # Ajdust ring radius by dx_km.
                        r0 = self.rho_km_vals[crange]
                        r = self.rho_km_vals[center]

                        # Compute psi for with stationary phase value
                        x = r-r0
                        x2 = np.square(x)
                        x3 = x2*x
                        x4 = np.square(x2)
                        x5 = x4*x
                        x6 = np.square(x3)

                        w_func = fw(x, w_init)
                    else:
                        crange += 1

                    z = x/d[center]
                    z2 = x2/d2[center]
                    z3 = x3/d3[center]
                    z4 = x4/d4[center]
                    z5 = x5/d5[center]
                    z6 = x6/d6[center]

                    psi_vals = b_0[center]-A_2[crange]*C_0[center]
                    psi_vals += z*(b_1[center]-A_2[crange]*C_1[center])
                    psi_vals += z2*(b_2[center]-A_2[crange]*C_2[center])
                    psi_vals += z3*(b_3[center]-A_2[crange]*C_3[center])
                    psi_vals += z4*(b_4[center]-A_2[crange]*C_4[center])
                    psi_vals += z5*(b_5[center]-A_2[crange]*C_5[center])
                    psi_vals += z6*(b_6[center]-A_2[crange]*C_6[center])

                    psi_vals *= z2*kD_vals[center]

                    # Compute kernel function for Fresnel inverse
                    if fwd:
                        ker = w_func*np.exp(1j*psi_vals)
                    else:
                        ker = w_func*np.exp(-1j*psi_vals)

                    # Range of diffracted data that falls inside the window
                    T = T_in[crange]

                    # Compute 'approximate' Fresnel Inversion for current point
                    T_out[center] = np.sum(ker*T)*self.dx_km*(1.0+1.0j)/(2.0*F)

                    # If normalization has been set, normalize the reconstruction
                    if self.norm:
                        T_out[center] *= __norm(self.dx_km, ker, F)
                    if self.verbose:
                        print(mes % (i, n_used, nw, loop), end="\r")
        if self.verbose:
            print("\n", end="\r")
        return T_out
