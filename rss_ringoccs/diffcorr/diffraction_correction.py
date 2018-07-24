"""
    Module Name:
        diffcorr
    Purpose:
        Provide functions and classes that aid in the process of
        Diffraction Correction / Fresnel Inversion. Additional
        functions for the purpose of forward modelling of
        reconstructed data and diffraction modelling are included.
        Several low-level functions that perform error checks for
        the main functions also exist, as well as functions that
        enable running shell scripts in Python.
"""

# Import dependencies for the diffcorr module
import time,os,platform
import numpy as np
from scipy.special import lambertw, iv
from scipy.constants import speed_of_light
from rss_ringoccs.tools import check_boole,check_real,check_pos_real

# Dictionary containing regions of interest within the Saturnian Rings.
region_dict = {
    'all'               : [65000.0,145000.0],
    'cringripples'      : [77690.0,77760.0],
    'encke'             : [132900.0,134200.0],
    'enckegap'          : [132900.0,134200.0],
    'janus'             : [96200.0,96800.0],
    'janusepimetheus'   : [96200.0,96800.0],
    'maxwell'           : [87410.0,87610.0],
    'maxwellringlet'    : [87410.0,87610.0],
    'titan'             : [77870.0,77930.0],
    'titanringlet'      : [77870.0,77930.0],
    'huygens'           : [117650.0,117950.0],
    'huygensringlet'    : [117650.0,117950.0]
    }

class rec_data(object):
    """
        Class:
            rec_data
        Purpose:
            Create an instance of a class that contains all of
            the input variables neccessary for diffraction
            reconstruction, including geometry, diffracted data,
            and calibration data.
        Arguments:
            NormDiff:
                An instance of the NormDiff class. This is a class
                that contains 8 variables necessary for diffraction
                correction:
                    rho_km_vals         - Ring Radius
                    p_norm_vals         - Diffracted Power
                    phase_rad_vals      - Diffracted Phase
                    B_rad_vals          - Ring Elevation
                    D_km_vals           - Spacecraft-RIP Distance
                    f_sky_hz_vals       - Sky Frequency
                    phi_rad_vals        - Ring Azimuth
                    rho_dot_kms_vals    - RIP Velocity
            res:
                Requested resolution for diffraction reconstruction.
                This is a floating point number that should be in 
                kilometers. Ex: res = 1.0 (1 km resolution).
            wtype:
                The name of the window function being applied during
                reconstruction. This is a string containing one of
                various allowed inputs. Ex: wtype = 'kb25'
        Keywords:
            bfac:
                A boolean for determining whether or not the 'b'
                factor will be used in the computation of window
                width. b is a function of rho_dot_kms_vals,
                sky_frequency, and the Allen deviation (see Sigma).
                Default is True. the b factor only play an important
                role when rho_dot_kms_vals is small. For zero Allen
                deviation, set bfac=False. Do NOT set sigma=0.
            sigma:
                A floating point number used to set what the Allen
                deviation is. If no value is set, the default value
                from the CASSINI spacecraft will be used (2e-13). 
                Spacecrafts with different Allen deviation need to
                use different values. For zero Allen deviation, set
                bfac=False. Do NOT set sigma=0.
        Attributes:
            res:
                Requested resolution for diffraction reconstruction.
                This is a floating point number that should be in 
                kilometers. Ex: res = 1.0 (1 km resolution).
            wtype:
                The name of the window function being applied during
                reconstruction. This is a string containing one of
                various allowed inputs. Ex: wtype = 'kb25'
            rho_km_vals:
                Ring radius of the ring intercept point. This is the
                intersection of the line from the spacecraft to Earth
                with the ring plane of the planet of interest. This
                is usually a numpy array containing floating
                point numbers. Values are in kilometers and are also
                nonnegative.
            p_norm_vals:
                Normalized power from the diffracted signal. This is
                a unitless variable from the normalization, and is
                usually an array of floating point numbers. Values
                are nonnegative.
            phase_rad_vals:
                The diffracted phase from the recieved signal. This
                is usually a numpy array containing floating point
                numbers. Values are between -2*Pi and 2*Pi and are
                in radians.
            B_rad_vals:
                The Ring Elevation Angle. This is the angle made
                between the line from the spacecraft to Earth and
                the ring plane of the planet of interest. This is
                usually a numpy array containing floating point
                numbers. Values range from -Pi to Pi and are in
                radians (Saturn values between -Pi/6 and Pi/6).
            D_km_vals:
                Distance from Spacecraft to Ring Intercept Point.
                The ring intercept point is the intersection of the
                line from the spacecraft to Earth and the ring plane
                of the planet of interest. D_km_vals is the distance,
                in kilometers, from the spacecraft to this point.
                This is usually a numpy array containing floating
                point numbers. Values are nonnegative.
            f_sky_hz_vals:
                Sky frequency. This is the frequency of the recieved
                signal, in Hertz. This is usually a numpy array
                containing floating point numbers. Value are in Hertz
                and are nonnegative.
            phi_rad_vals:
                The observed ring azimuth angle. This is the angle
                made by the ring intercept point and the x-axis of
                the ring plane. The x-axis is defined as the line
                parallel to the projection of the spacecraft-Earth
                line onto the ring plane, which contains the center
                of the planet of interest. This is usually a numpy
                array containing floating point numbers. Values are
                in radians are range from -2*Pi to 2*Pi.
            rho_dot_kms_vals:
                The time derivative of the ring intercept point. The
                ring intercept point is the intersection of the line
                from the spacecraft to Earth with the ring plane of
                the planet of interest. This is usually a numpy array
                containing floating point numbers. For proper
                reconstruction, these values must be positive. Thus,
                during observation with negative rho_dot_kms_vals,
                one must swap all variables in time so that the
                rho_dot_kms_vals variable is positive. Observations
                of positive and negative rho_dot_kms_vals must be
                split. Reconstruction is not possible when
                rho_dot_kms_vals is zero. Values are in km/s.
            T_hat_vals:
                Complex transmittance of the diffracted signal. This
                is a numpy array containing floating point complex
                numbers. Defined as T_hat = sqrt(Power)*Exp(i*Phase)
            F_km_vals:
                The Fresnel scale. This is a length that is intrinsic
                to the reconstruction of the diffracted signal. For
                more information, please see the references below.
                This is usually a numpy array containing floating
                point numbers. Value are in kilometers and are
                nonnegative.
            w_km_vals:
                Window width as a function of ring intercept point.
                This is the width of the window function needed to
                compute reconstruction at the desired resolution.
                This is usually a numpy array of floating point
                numbers. Values are in kilometers are nonnegative.
            mu_vals:
                The sine of the absolute value of the B_rad_vals
                attribute. Numpy array that is unitless.
            lambda_sky_km_vals:
                The wavelength that corresponds to the f_sky_Hz_vals
                attribute. Numpy array that is in kilometers and
                nonnegative.
            dx_km:
                The sampling spacing. This is the distance between
                points in the ring plane that were sampled by the
                observation. This is a single floating point number.
            norm_eq:
                Normalized equivalent width. A certain mathematical
                value used in the evaluation of window widths and to
                represent equivalent widths of windows. This is a
                single floating point number that is unitless.
            t_oet_spm_vals:
                Observed Event Time. The time the signal was recorded
                on Earth. This is usually a numpy array containing
                floating point numbers. Values are in seconds.
                Ring Event Time. The time the signal was diffracted
                by the rings. This is usually a numpy array
                containing floating point numbers. Values are in
                seconds.
            t_set_spm_vals:
                Spacecraft Event Time. The time the singal was
                emitted from the spacecraft. This is usually a numpy
                array containing floating point numbers. Values are
                in seconds.
            rho_corr_pole_km_vals:
                Ring radius corrections based on the planet's
                pole direction. This is usually a numpy array
                containing floating point numbers. Values are in
                kilometers.
            rho_corr_timing_km_vals:
                Ring radius correction from timing offsets. This is
                usually a numpy array containing floating point
                numbers. Values are in kilometers.
            phi_rl_rad_vals:
                Observed Ring Longitude. A numpy array containing
                floating point numbers that vary between -2*Pi and
                2*Pi.
        Dependencies:
            [1] sys
            [2] numpy
            [3] rss_ringoccs.tools subpackage
            [4] scipy.special
            [5] scipy.constants
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] S. W. Asmar, R. G. French, E. A. Marouf, P. Schinder,
                J. W. Armstrong, P. Tortora, L. Iess, A. Anabtawi,
                A. J. Kliore, M. Parisi, M. Zannoni, and D. Kahan,
                Cassini Radio Science User's Guide, September, 2014,
                https://pds-rings.seti.org/cassini/rss/
    
	"""

    def __init__(self, NormDiff,res,wtype,bfac=True,sigma=2.e-13):

        # Perform an error check to ensure inputs are valid
        if not check_pos_real(res):
            raise TypeError("res must be a positive number")
        if not (type(wtype)==type("Hi!")):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        if not check_boole(bfac):
            raise TypeError("bfac must be Boolean: True/False")
        if (not check_pos_real(sigma)):
                raise TypeError("sigma must be a positive number.")

        # Define/list all attributes of the class, setting to None.
        self.res                        = None  # Resolution (km)
        self.wtype                      = None  # Window Type
        self.rho_km_vals                = None  # Ring Radius (km)
        self.p_norm_vals                = None  # Normalized Power
        self.phase_rad_vals             = None  # Phase (Radians)
        self.B_rad_vals                 = None  # Elevation Angle (Radians)
        self.D_km_vals                  = None  # RIP Distance (km)
        self.f_sky_hz_vals              = None  # Signal Frequency (Hz)
        self.phi_rad_vals               = None  # Azimuthal Angle (Radians)
        self.rho_dot_kms_vals           = None  # RIP Velocity (km/s)
        self.T_hat_vals                 = None  # Complex Transmittance
        self.F_km_vals                  = None  # Fresnel Scale (km)
        self.w_km_vals                  = None  # Window Width (km)
        self.mu_vals                    = None  # sin(|B|)
        self.lambda_sky_km_vals         = None  # Signal Wavelength (km)
        self.dx_km                      = None  # Sampling Distance (km)
        self.norm_eq                    = None  # Normalized Equivalent Width
        self.t_oet_spm_vals             = None  # Observed Event Time (seconds)
        self.t_ret_spm_vals             = None  # Ring Event Time (seconds)
        self.t_set_spm_vals             = None  # Spacecraft ET (seconds)
        self.rho_corr_pole_km_vals      = None  # RIP Pole Corrections (km)
        self.rho_corr_timing_km_vals    = None  # RIP Timing corrections (km)
        self.phi_rl_rad_vals            = None  # Ring Longitude (Radians)

        # Retrieve variables from the NormDiff class, setting as attributes.
        self.res                        = res
        self.wtype                      = wtype.replace(" ", "").lower()
        self.rho_km_vals                = NormDiff.rho_km_vals
        self.p_norm_vals                = NormDiff.p_norm_vals
        self.phase_rad_vals             = -NormDiff.phase_rad_vals
        self.B_rad_vals		            = NormDiff.B_rad_vals
        self.D_km_vals		            = NormDiff.D_km_vals
        self.f_sky_hz_vals              = NormDiff.f_sky_hz_vals
        self.phi_rad_vals               = NormDiff.phi_rad_vals
        self.rho_dot_kms_vals           = NormDiff.rho_dot_kms_vals
        self.t_ret_spm_vals             = NormDiff.t_ret_spm_vals         
        self.t_set_spm_vals             = NormDiff.t_set_spm_vals         
        self.rho_corr_pole_km_vals      = NormDiff.rho_corr_pole_km_vals
        self.rho_corr_timing_km_vals    = NormDiff.rho_corr_timing_km_vals
        self.phi_rl_rad_vals            = NormDiff.phi_rl_rad_vals

        # Compute mu: Sin(|B|)
        self.mu_vals = np.sin(np.abs(self.B_rad_vals))

        # Try the two names attributed to Observed Event Time.
        try: self.t_oet_spm_vals = NormDiff.t_oet_spm_vals         
        except AttributeError:
            try: self.t_oet_spm_vals = NormDiff.spm_vals
            except AttributeError: pass

        # Check that rho_km_vals is increasing and the rev isn't a chord occ.
        drho = [np.min(self.rho_dot_kms_vals),np.max(self.rho_dot_kms_vals)]
        dx   = self.rho_km_vals[1]-self.rho_km_vals[0]
        if (drho[0] < 0) and (drho[1] > 0):
            raise ValueError("drho/dt had positive and negative values.")
        if (dx > 0) and (drho[1] < 0):
            self.rho_dot_kms_vals=np.abs(self.rho_dot_kms_vals)
        if (dx < 0):
            self.rho_km_vals      = self.rho_km_vals[::-1]
            self.phase_rad_vals   = self.phase_rad_vals[::-1]
            self.p_norm_vals      = self.p_norm_vals[::-1]
            self.phi_rad_vals     = self.phi_rad_vals[::-1]
            self.B_rad_vals       = self.B_rad_vals[::-1]
            self.f_sky_hz_vals    = self.f_sky_hz_vals[::-1]
            self.D_km_vals        = self.D_km_vals[::-1]
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])

        # Check that all variables from NormDiff are the same size.
        n_rho = np.size(self.rho_km_vals)
        if (np.size(self.phase_rad_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(phase) != len(rho)")
        if (np.size(self.p_norm_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(power) != len(rho)")
        if (np.size(self.phi_rad_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(phi) != len(rho)")
        if (np.size(self.B_rad_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(B) != (rho)")
        if (np.size(self.f_sky_hz_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(frequency) != len(rho)")
        if (np.size(self.D_km_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(D) != len(rho)")
        if (np.size(self.rho_dot_kms_vals) != n_rho):
            raise ValueError("Bad NormDiff: len(rho_dot) != len(rho)")

        # Compute necessary variables for diffraction correction.
        self.lambda_sky_km_vals = speed_of_light*0.001 / self.f_sky_hz_vals
        self.dx_km              = self.rho_km_vals[1] - self.rho_km_vals[0]
        self.T_hat_vals         = np.sqrt(np.abs(self.p_norm_vals))*\
            np.exp(1j*self.phase_rad_vals)

        # Compute the Fresnel Scale (See MTR86 Equation 6)
        cb                      = np.cos(self.B_rad_vals)
        sb                      = np.sin(self.B_rad_vals)
        sp                      = np.sin(self.phi_rad_vals)
        self.F_km_vals          = np.sqrt(0.5 * self.lambda_sky_km_vals *\
            self.D_km_vals * (1 - (cb*cb) * (sp*sp)) / (sb*sb))
        
        # Compute the Normalized Equaivalent Width (See MTR86 Equation 20)
        self.norm_eq            = self.__func_dict[wtype]

        # Compute the window width. (See MTR86 Equations 19, 32, and 33).
        if bfac:
            omega = 2.0 * np.pi * self.f_sky_hz_vals
            alpha = (omega*omega) * (sigma*sigma)/(2.0 * self.rho_dot_kms_vals)
            P     = res / (alpha * (self.F_km_vals*self.F_km_vals))
            # The inverse exists only if P>1.
            if (np.min(P) < 1.0001):
                raise ValueError("\nWarning: Bad Points!\n\
                \rEither rho_dot_kms_vals, F_km_vals, or res_km is too\n\
                \rsmall. Exclude points or request coarser resolution.")
            self.w_km_vals = self.norm_eq*np.abs(((P-1)*\
                lambertw(np.exp(P/(1-P))*P/(1-P))+P)/(P-1)/alpha)
        else:self.w_km_vals=2.0*self.norm_eq*self.F_km_vals*self.F_km_vals/res
    
    # Dictionary containing pre-computed normalized equivalent widths.
    __func_dict = {
        "rect"      : 1.00000000,
        "coss"      : 1.50000000,
        "kb20"      : 1.49634231,
        "kb25"      : 1.65191895,
        "kb35"      : 1.92844639,
        "kbmd20"    : 1.52048174,
        "kbmd25"    : 1.65994218}

class diffraction_correction(object):

    def __init__(self,dat,res,rng="all",wtype="kb25",fwd=False,
        norm=True,verbose=False,bfac=True,fft=False,psitype="full"):
        t1       = time.time()

        if not check_pos_real(res):
            raise TypeError("res must be a positive number")
        if not (type(wtype)==type("Hi!")):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        if not check_boole(fwd):
            raise TypeError("fwd must be Boolean: True/False")
        if not check_boole(norm):
            raise TypeError("norm must be Boolean: True/False")
        if not check_boole(bfac):
            raise TypeError("bfac must be Boolean: True/False")
        if not check_boole(fft):
            raise TypeError("fft must be Boolean: True/False")
        if not check_boole(verbose):
            raise TypeError("verbose must be Boolean: True/False")
        if not (type(psitype)==type("Hi!")):
            raise TypeError("psitype must be a string. Ex: 'full'")
        tr = type(rng)
        if (not check_real(rng)) and (tr != type([1])) and (tr != type("Hi")):
            raise TypeError("rng must be of the form [a,b] or a valid string.")
        else: del tr

        # Define/list all attributes of the class, setting to None.
        self.res                = None  # Resolution (km)
        self.wtype              = None  # Window Type
        self.rng                = None  # Reconstructed Range (km)
        self.rngreq             = None  # Requested Range
        self.rho_km_vals        = None  # Ring Radius (km)
        self.p_norm_vals        = None  # Normalized Power
        self.phase_rad_vals     = None  # Phase (Radians)
        self.B_rad_vals         = None  # Elevation Angle (Radians)
        self.D_km_vals          = None  # RIP Distance (km)
        self.f_sky_hz_vals      = None  # Signal Frequency (Hz)
        self.phi_rad_vals       = None  # Azimuthal Angle (Radians)
        self.rho_dot_kms_vals   = None  # RIP Velocity (km/s)
        self.T_hat_vals         = None  # Complex Transmittance
        self.F_km_vals          = None  # Fresnel Scale (km)
        self.w_km_vals          = None  # Window Width (km)
        self.mu_vals            = None  # sin(|B|)
        self.lambda_sky_km_vals = None  # Signal Wavelength (km)
        self.dx_km              = None  # Sampling Distance (km)
        self.norm_eq            = None  # Normalized Equivalent Width
        self.n_used             = None  # Number of points used
        self.start              = None  # Starting point
        self.finish             = None  # Final point
        self.T_vals             = None  # Reconstructed Transmittance
        self.power_vals         = None  # Reconstructed power
        self.tau_vals           = None  # Reconstructed optical depth
        self.phase_vals         = None  # Reconstructed phase
        self.p_norm_fwd_vals    = None  # Forward Power
        self.T_hat_fwd_vals     = None  # Forward Transmittance
        self.phase_fwd_vals     = None  # Forward Phase
        self.verbose            = None  # Use of verbose (Boolean)
        self.norm               = None  # Normalization (Boolean)
        self.fwd                = None  # Forward (Boolean)
        self.fft                = None  # Use of FFT's (Boolean)
        self.bfac               = None  # Use of b factor (Boolean)
        self.psitype            = None  # Psitype used (String)
        self.history            = None  # History of Inputs
        self.dathist            = None  # History from rec_data

        # Assign keywords/arguments to their corresponding attribute.
        self.res        = res           # Resolution (km)
        self.wtype      = wtype         # Window Type (String)
        self.rngreq     = rng           # Requested Range
        self.norm       = norm          # Normalization (Boolean)
        self.verbose    = verbose       # Use of verbose (Boolean)
        self.fwd        = fwd           # Forward (Boolean)
        self.fft        = fft           # Use of FFT's (Boolean)
        self.bfac       = bfac          # Use of b factor (Boolean)
        self.psitype    = psitype       # Psitype used (String)
        self.dathist    = dat.history   # rec_data History

        # Create an instance of the rec_data class, containing data.
        recdata = rec_data(dat, res, wtype, bfac=bfac)

        # Save the attributes from the rec_data instance.
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.B_rad_vals         = recdata.B_rad_vals
        self.D_km_vals          = recdata.D_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq

        # Remove rec_data instance and keywords for memory and clarity.
        del recdata

        # From the requested range, extract array of the form [a,b]
        if (type(rng) == type('Hello')):
            rng = rng.replace(" ","").lower()
            if (rng in region_dict): self.rng = np.array(region_dict[rng])
            else:
                erm = ""
                for key in region_dict: erm = "%s%s\n" % (erm,key)
                raise ValueError("Illegal Range. Allowed Strings:\n%s"%erm)
        elif (np.size(rng) < 2):
            raise TypeError("Only one value given for range. Use rng = [a,b]")
        else:   self.rng = np.array([np.min(rng),np.max(rng)])

        # Compute the starting point and the number of points used.
        rho         = self.rho_km_vals          # Ring radius.
        w_max       = np.max(self.w_km_vals)    # Largest window used.

        # Compute the smallest and largest allowed radii for reconstruction.
        rho_min_lim = np.min(rho)+np.ceil(w_max/2.0)
        rho_max_lim = np.max(rho)-np.ceil(w_max/2.0)

        # Compute the smallest and largest values within requested range.
        rho_start   = rho[np.min((rho >= np.min(self.rng)).nonzero())]
        rho_end     = rho[np.max((rho <= np.max(self.rng)).nonzero())]

        # Compute the start and end point for reconstruction.
        rho_min     = np.max([rho_min_lim,rho_start])
        rho_max     = np.min([rho_max_lim,rho_end])
        self.start  = int(np.min((rho >= rho_min).nonzero()))
        self.finish = int(np.max((rho <= rho_max).nonzero()))
        self.n_used = 1 + (self.finish - self.start)

        # Delete unnecessary variables for clarity and memory.
        del rho,w_max,rho_min_lim,rho_max_lim,rho_start,rho_end
        del rho_min,rho_max,rng,norm,fwd,fft,bfac,psitype,verbose

        # self.__trim_inputs()
        self.T_vals = self.__finv()
        if self.fwd:
            self.T_hat_fwd_vals     = self.__ffwd()
            self.p_norm_fwd_vals    = power_func(self.T_hat_fwd_vals)
            self.phase_fwd_vals     = phase_func(self.T_hat_fwd_vals)

        self.power_vals = np.abs(self.T_vals*self.T_vals)
        self.phase_vals = -np.arctan2(np.imag(self.T_vals),np.real(self.T_vals))
        crange          = (self.power_vals>0).nonzero()
        tau             = np.zeros(np.size(self.power_vals))
        tau[crange]     = -self.mu_vals[crange]*\
            np.log(np.abs(self.power_vals[crange]))
        self.tau_vals   = tau

        self.__trim_attributes(self.fwd)

        self.history = self.__write_hist_dict()

        t2 = time.time()
        print("Computation Time: ",t2-t1,end="\r")

    # Definition of the Rectangular window function.
    def __rect(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute window function
        w_func = np.zeros(nw_pts) + 1.0
        return w_func

    # Definition of the Squared Cosine window function.
    def __coss(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Compute window function.
        w_func = np.cos(np.pi * x / w_in)**2
        return w_func

    # Definition of the Kaiser-Bessel 2.0 window function.
    def __kb20(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb20 is 2.0
        alpha  = 2.0 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    # Definition of the Kaiser-Bessel 2.5 window function.
    def __kb25(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb25 is 2.5.
        alpha  = 2.5 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    # Definition of the Kaiser-Bessel 3.5 window function.
    def __kb35(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb35 is 3.5.
        alpha  = 3.5 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    # Modified Kaiser-Bessel 2.0 window function.
    def __kbmd20(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kbmd20 is 2.0.
        alpha  = 2.0*np.pi
        # Compute window function.
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    # Modified Kaiser-Bessel 2.5 window function.
    def __kbmd25(w_in, dx):
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kbmd25 is 2.5.
        alpha  = 2.5*np.pi
        # Compute window function.
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    # Function dictionary with normalized equivalent widths.
    __func_dict = {
        "__rect" : {"func" : __rect, "normeq" : 1.00000000},
        "__coss" : {"func" : __coss, "normeq" : 1.50000000},
        "__kb20" : {"func" : __kb20, "normeq" : 1.49634231},
        "__kb25" : {"func" : __kb25, "normeq" : 1.65191895},
        "__kb35" : {"func" : __kb35, "normeq" : 1.92844639},
        "__kbmd20" : {"func" : __kbmd20, "normeq" : 1.52048174},
        "__kbmd25" : {"func" : __kbmd25, "normeq" : 1.65994218}
        }

    def __trim_attributes(self,fwd):
        # Get rid of uncomputed values and keep only what was processed.
        start  = self.start                 # Starting Point.
        n_used = self.n_used                # Number of point used.
        crange = np.arange(n_used)+start    # Processed range.

        self.rho_km_vals         = self.rho_km_vals[crange]
        self.p_norm_vals         = self.p_norm_vals[crange]
        self.phase_rad_vals      = self.phase_rad_vals[crange]
        self.B_rad_vals          = self.B_rad_vals[crange]
        self.D_km_vals           = self.D_km_vals[crange]
        self.f_sky_hz_vals       = self.f_sky_hz_vals[crange]
        self.phi_rad_vals        = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals    = self.rho_dot_kms_vals[crange]
        self.T_hat_vals          = self.T_hat_vals[crange]
        self.F_km_vals           = self.F_km_vals[crange]
        self.w_km_vals           = self.w_km_vals[crange]
        self.mu_vals             = self.mu_vals[crange]
        self.lambda_sky_km_vals  = self.lambda_sky_km_vals[crange]
        self.T_vals              = self.T_vals[crange]
        self.power_vals          = self.power_vals[crange]
        self.tau_vals            = self.tau_vals[crange]
        self.phase_vals          = self.phase_vals[crange]

        #If the forward model was run, trim those attributes as well.
        if fwd:
            self.p_norm_fwd_vals = self.p_norm_fwd_vals[crange]
            self.T_hat_fwd_vals  = self.T_hat_fwd_vals[crange]
            self.phase_fwd_vals  = self.phase_fwd_vals[crange]
    
    def __trim_inputs(self):
        start   = self.start
        n_used  = self.n_used
        rho     = self.rho_km_vals
        rstart  = rho[start]
        rfin    = rho[start+n_used+1]
        w_vals  = self.w_km_vals
        w       = np.ceil(np.max(w_vals[start:start+n_used+1])/2.0)
        nst = np.min((rho>=(rstart-w)).nonzero())
        nen = np.max((rho<=(rfin+w)).nonzero())

        del n_used, rho, rstart, rfin, w_vals, w

        nreq   = 1 + (nen - nst)
        crange = np.arange(nst)+nreq
        self.rho_km_vals         = self.rho_km_vals[crange]
        self.start               = start-nreq
        self.p_norm_vals         = self.p_norm_vals[crange]
        self.phase_rad_vals      = self.phase_rad_vals[crange]
        self.B_rad_vals          = self.B_rad_vals[crange]
        self.D_km_vals           = self.D_km_vals[crange]
        self.f_sky_hz_vals       = self.f_sky_hz_vals[crange]
        self.phi_rad_vals        = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals    = self.rho_dot_kms_vals[crange]
        self.T_hat_vals          = self.T_hat_vals[crange]
        self.F_km_vals           = self.F_km_vals[crange]
        self.w_km_vals           = self.w_km_vals[crange]
        self.mu_vals             = self.mu_vals[crange]
        self.lambda_sky_km_vals  = self.lambda_sky_km_vals[crange]

        del nreq, crange

    # Definition of the Approximate Fresnel Inverse (MTR86 Equation 15)
    def __fresinv(self,T_hat,ker,dx,f_scale):
        T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
        return T

    # FFT Approximation of Fresnel Inverse (MTR86 Equation 15)
    def __fresinvfft(self,T_hat,ker,dx,f_scale):
        nw = np.size(T_hat)
        fft_t_hat       = np.fft.fft(T_hat)
        fft_conv        = np.fft.fft(ker)
        inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
        inv_t_hat      *= dx*(np.complex(1.0,1.0))/(2.0*f_scale)
        T               = inv_t_hat[int((nw-1)/2)]
        return T

    # Window Normalization Function.
    def __normalize(self,dx,w_func,psi_vals,f_scale):
        ker       = np.exp(-1j * psi_vals)              # Kernel Function
        T1        = np.abs(np.sum(w_func * ker) * dx)   # Freespace Integral
        norm_fact = np.sqrt(2.0) * f_scale / T1         # Normalization Factor
        return norm_fact
    
    # 'Do Nothing' Function for when norm=False is set.
    def __normalize_do_nothing(self,w,x,y,z):
        return 1

    # Printing function for when verbose=True is set.
    def __verbose_print(self,i,n_used,nw,loop):
        print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
            % (i,n_used,nw,loop),end="\r")
    
    # 'Do Nothing' function for when verbose=False is set.
    def __verbose_do_nothing(self,i,n_used,nw,loop):
        pass

    def __finv(self):
        # Retrieve variables.
        w_vals          = self.w_km_vals
        rho_vals        = self.rho_km_vals
        phi_rad_vals    = self.phi_rad_vals
        d_vals          = self.D_km_vals
        B_rad_vals      = self.B_rad_vals
        lambda_vals     = self.lambda_sky_km_vals
        T_hat_vals      = self.T_hat_vals
        F_vals          = self.F_km_vals
        fft             = self.fft
        verbose         = self.verbose
        wtype           = "%s%s" % ("__",self.wtype)
        start           = self.start
        n_used          = self.n_used
        dx              = self.dx_km
        # Compute necessary variables.
        kD_vals         = 2.0 * np.pi * d_vals / lambda_vals
        cosb            = np.cos(B_rad_vals)
        cosb_D          = cosb/d_vals
        cosb2           = cosb*cosb
        cosphi0         = np.cos(phi_rad_vals)
        sinphi0         = np.sin(phi_rad_vals)
        dsq             = d_vals*d_vals
        rsq             = rho_vals*rho_vals
        # Define functions
        fw              = self.__func_dict[wtype]["func"]
        if fft:         finv    = self.__fresinvfft
        else:           finv    = self.__fresinv
        if self.norm:   nrm     = self.__normalize
        else:           nrm     = self.__normalize_do_nothing
        if verbose:     vrb     = self.__verbose_print
        else:           vrb     = self.__verbose_do_nothing
        # Calculate the corrected complex amplitude, point by point
        T_vals     = T_hat_vals * 0.0
        center     = start
        w_init     = w_vals[center]
        w_func     = fw(w_init,dx)
        nw         = np.size(w_func)
        crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))-1
        phi_s_rad  = phi_rad_vals[center]
        cp         = np.cos(phi_s_rad)
        sp         = np.sin(phi_s_rad)
        dphi_fac   = ((cosb2)*cosphi0*sinphi0/(1.0-(cosb2)*(sinphi0*sinphi0)))
        for i in np.arange(n_used):
            r0          = rho_vals[center]
            r02         = rsq[center]
            d2          = dsq[center]
            cb_d        = cosb_D[center]
            cb2         = cosb2[center]
            cp0         = cosphi0[center]
            sp0         = sinphi0[center]
            kD          = kD_vals[center]
            w           = w_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(center-(nw-1)/2,1+center+(nw-1)/2)
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = dphi_fac[center] * (r - r0) / r0
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange    +=1
                r          = rho_vals[crange]
                r2         = rsq[crange]
                xi         = cb_d * (r0 * cp0 - r * cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = cb_d * r * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = (2.0*v1 + v2) / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            loop = 0
            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
                xi         = cb_d * (r0 * cp0 - r * cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = cb_d * r * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = (2.0*v1 + v2) / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
 
            # Compute psi and then compute the forward model.
            xi       = cb_d * (r0 * cp0 - r * cp)
            eta      = ((r02) + (r2) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / d2
            psi_vals = kD * (np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi))
            F        = F_vals[center]
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center]   = finv(T_hat,ker,dx,F)
            T_vals[center]  *= nrm(dx,w_func,psi_vals,F)
            center          += 1
            vrb(i,n_used,nw,loop)
        self.phi_s = phi_s_rad
        return T_vals

    def __write_hist_dict(self):
        """
        This creates a history dictionary.

        Returns:
            geo_hist (dict): Dictionary with "user name", "host name",
                    "run date", "python version", "operating system",
                    "source file", "input variables", and "input keywords".
        """
        user_name = os.getlogin()
        host_name = os.uname()[1]
        run_date  = time.ctime() + ' ' + time.tzname[0]
        python_v  = platform.python_version()
        opsys     = os.uname()[0]
        src_file  = __file__.split('/')[-1]
        src_dir   = __file__.rsplit('/',1)[0] +'/'
        rngreq    = (np.min(self.rng),np.max(self.rng))

        tau_hist = {
            "User Name"         : user_name,
            "Host Name"         : host_name,
            "Run Date"          : run_date,
            "Python Version"    : python_v,
            "Operating System"  : opsys,
            "Source Directory"  : src_dir,
            "Source File"       : src_file,
            "Input Variables"   : {
                "NormDiff Dictionary": self.dathist},
            "Input Keywords"    : {
                "b Factor"                      : self.bfac,
                "Resolution (km)"               : self.res,
                "Tapering Function"             : self.wtype,
                "Requested Radius (km)"         : rngreq,
                "Normalization"                 : self.norm,
                "Forward Model Data"            : self.fwd,
                "FFT's Used"                    : self.fft,
                "Psi Approximation"             : self.psitype},
        }
        return tau_hist