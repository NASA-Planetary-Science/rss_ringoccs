# Import dependencies for the diffcorr module
import time
import sys
import numpy as np
from scipy.special import lambertw, iv
from rss_ringoccs.tools.write_history_dict import write_history_dict
SPEED_OF_LIGHT_KM = 299792.4580

# Dictionary containing regions of interest within the Saturnian Rings.
region_dict = {
    'all':             [65000.0, 145000.0],
    'cringripples':    [77690.0, 77760.0],
    'encke':           [132900.0, 134200.0],
    'enckegap':        [132900.0, 134200.0],
    'janusepimetheus': [96200.0, 96800.0],
    'maxwell':         [87410.0, 87610.0],
    'maxwellringlet':  [87410.0, 87610.0],
    'titan':           [77870.0, 77930.0],
    'titanringlet':    [77870.0, 77930.0],
    'huygens':         [117650.0, 117950.0],
    'huygensringlet':  [117650.0, 117950.0]
    }


class DiffractionCorrection(object):
    """
        Class:
            diffraction_correction
        Purpose:
            Perform diffraction correction for a ring occultation
            on a data set that is a near radially symmetric function
            of the ring radius, or ring intercept point (RIP).
        Arguments:
            dat:    
                The data set, usually an instance of the NormDiff
                class from the rss_ringoccs Calibration subpackage.
                This instance MUST contain the following attributes
                and MUST have the same names.
                    rho_km_vals:      Ring Radius (km)
                    phi_rad_vals:     Ring Azimuth Angle (Radians)
                    p_norm_vals:      Normalized Power
                    phase_rad_vals:   Phase (Radians)
                    B_rad_vals:       Elevation Angle (Radians)
                    D_km_vals:        RIP-Distance (km)
                    f_sky_hz_vals:    Sky Frequency (Hertz)
                    rho_dot_kms_vals: RIP-velocity (km/s)
                    history:          History dictionary
            res:    
                The requested resolution for processing (km). This
                must be a positive real number, that is, a positive
                floating point number or integer.
        Keywords:
            rng:    
                The request range for diffraction correction.
                Preferred input is rng = [a,b]. Arrays are
                allowed and the range will be set as:
                    rng = [MIN(array),MAX(array)]
                Finally, certain strings containing a few of the
                regions of interests within the rings of Saturn
                are allowed. Permissible strings are:
                    'maxwell', 'titan', 'huygens', and 'encke'.
                Strings are neither case nor space sensitive.
                For other planets use rng = [a,b]. Default value
                is set to 'all' which processing [65,000,140,000]
                Values MUST be set in kilometers.
            wtype:  
                The requested tapering function for diffraction
                correction. A string with several allowed inputs:
                    'rect'      Rectangular Window.
                    'coss'      Squares Cosine Window.
                    'kb20'      Kaiser-Bessel 2.0 Window.
                    'kb25'      Kaiser-Bessel 2.5 Window.
                    'kb35'      Kaiser-Bessel 3.5 Window.
                    'kbmd20'    Modified kb20 Window.
                    'kbmd25'    Modified kb25 Window.
                The variable is neither case nor space sensitive.
                Default window is set to 'kb25'
            fwd:    
                A Boolean for determining whether or not
                forward modelling will be computed. This is good
                starting point for deciding if the diffraction
                correction is physically significant or valid. If
                the reconstruction is good, the forward model
                should reproduce the p_norm_vals attribute from
                the input dat instance. Default is set to False.
            norm:
                A Boolean for determining whether or not the
                reconstructed complex transmittance is normalize
                by the window width. This normalization is the
                complex transmittance that is computed by using
                free space divided by the complex transmittance
                that is computed using free space weighted by the
                selected tapering function. Default is True.
            bfac:
                A Boolean for determining whether or not the
                'b' factor in the window width computation is
                used. This is equivalent to setting the Allen
                Deviation from the spacecraft to a positive value
                or to zero. If set to False, the Allen Deviation
                is assumed to be zero. If set to True the Allen
                Deviation is set to 2e-13. If set to a positive
                real number, the Allen Deviation will be assumed
                to be that real number. Default is True.
            fft:    
                A Boolean for determining whether or not FFT's will
                be used for computing the complex transmittance. The
                use of FFT's assumes that the geometry of the system
                is such that the integral that is used to compute the
                complex transmittance is of the form of a
                convolution, and that the convolution theorem may be
                applied to it. Default is set to False.
            psitype:
                A string for determining what approximation to the
                geometrical 'psi' function is used. Several strings
                are allowed:
                    'full'      No Approximation is applied.
                    'taylor2'   Second order Taylor Series.
                    'taylor3'   Third order Taylor Series.
                    'taylor4'   Fourth order Taylor Series.
                    'MTR2'      Second Order Series from MTR86.
                    'MTR3'      Third Order Series from MTR86.
                    'MTR4'      Fourth Order Series from MTR86.
                The variable is neither case nor space sensitive.
                Default is set to 'full'.
            verbose:
                A Boolean for determining if various pieces of
                information are printed to the screen or not.
                Default is False.
        Outputs:
            T_hat_vals:
                Complex transmittance of the diffracted data.
            F_km_vals:
                Fresnel scale (km).
            w_km_vals:
                Window width as a function of ring radius (km).
            mu_vals:
                The sine of the elevation angle.
            lambda_sky_km_vals:
                Wavelength of the recieved signal (km).
            dx_km:  
                Radial spacing between points (km).
            norm_eq:        
                Normalized equivalent width of the window function.
            n_used:
                Number of points that were processed (integer).
            start:
                Starting point used for processing (integer).
            T_vals:
                Complex transmittance of reconstructed data.
            power_vals:
                Normalized power of the reconstructed data.
            tau_vals:
                Normalized optical depth of the reconstructed data.
            phase_vals:
                Phase of the reconstructed data (Radians).
            p_norm_fwd_vals:
                Normalized power of forward model (fwd=True).
            T_hat_fwd_vals:
                Complex transmittance of forward model (fwd=True).
            phase_fwd_vals:
                Phase of forward model (fwd=True).
            history:
                History dictionary of the runtime OS settings.
        Dependencies:
            [1] numpy
            [2] scipy
            [3] diffcorr
            [4] 
        Notes:
            [1] 
        References:

        Examples:

        History:
            Created: RJM - 2018/05/16 5:40 P.M.
    """
    def __init__(self, NormDiff, res, rng="all", wtype="kb25", fwd=False,
                 norm=True, verbose=False, bfac=True, sigma=2.e-13,
                 fft=False, psitype="full"):
        t1 = time.time()

        if verbose:
            print("Processing Diffraction Correction:")
            print("\tRunning Error Check on Input Arguments...")
        if (not isinstance(res, float)) or (not isinstance(res, int)):
            try:
                res = float(res)
            except TypeError:
                raise TypeError("res must be a positive floating point number")
        if (res <= 0.0):
            raise ValueError("res must be a positive floating point number")
        if not isinstance(wtype, str):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        elif not (wtype in self.__func_dict):
            erm = ""
            for key in self.__func_dict:
                erm = "%s%s\n" % (erm, key)
            raise ValueError("\nIllegal string used for wtype.\n"
                             "\rAllowed Strings:\n%s" % erm)
        if not isinstance(fwd, bool):
            raise TypeError("fwd must be Boolean: True/False")
        if not isinstance(norm, bool):
            raise TypeError("norm must be Boolean: True/False")
        if not isinstance(bfac, bool):
            raise TypeError("bfac must be Boolean: True/False")
        if not isinstance(fft, bool):
            raise TypeError("fft must be Boolean: True/False")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be Boolean: True/False")
        if not isinstance(psitype, str):
            raise TypeError("psitype must be a string. Ex: 'full'")
        if (not isinstance(rng, str)) and (not isinstance(rng, list)):
            try:
                rng = [np.min(rng), np.max(rng)]
                if (np.min(rng) < 0):
                    raise ValueError(
                        "Minimum requested range must be positive"
                        )
            except TypeError:
                raise TypeError("rng must be a list [a,b] or a valid string")
        elif isinstance(rng, list):
            if (not all(isinstance(x, float) for x in rng)):
                try:
                    for i in range(np.size(rng)):
                        rng[i] = float(rng[i])
                except TypeError:
                    raise TypeError(
                        "rng must be a list of floating point numbers"
                        )
            elif (np.size(rng) > 2):
                raise TypeError("rng must contain two numbers: rng = [a,b]")
            elif (np.min(rng) < 0.0):
                raise ValueError("Minimum requested range must be positive")
        elif isinstance(rng, str):
            rng = rng.replace(" ", "").lower()
            if not (rng in region_dict):
                erm = ""
                for key in region_dict:
                    erm = "%s%s\n" % (erm, key)
                raise ValueError("\nIllegal string used for rng.\n"
                                 "\rAllowed Strings:\n%s" % erm)
        if (not isinstance(sigma,float)):
            try:
                sigma = float(sigma)
            except TypeError:
                raise TypeError(
                    "sigma must be a positive floating point number."
                    )
        if (np.min(sigma) < 0.0):
            raise TypeError("sigma must be a positive floating point number.")

        if verbose:
            print("\tAssigning Inputs as Attributes...")
        # Set forward power variable to None in case it isn't defined later.
        self.p_norm_fwd_vals = None

        # Set forward phase and power to None as well.
        self.T_hat_fwd_vals = None
        self.phase_fwd_vals = None

        # Assign resolution, forward, and FFT variables as attributes.
        self.res = res
        self.fwd = fwd
        self.fft = fft

        # Assing window type and Allen deviation variables as attributes.
        self.wtype = wtype.replace(" ", "").lower()
        self.sigma = sigma

        # Assign normalization and b-factor variables as attributes.
        self.norm = norm
        self.bfac = bfac

        # Assign requested range variable as an attribute.
        self.rngreq = rng

        # Assign verbose and psi approximation variables as attributes.
        self.verbose = verbose
        self.psitype = psitype.replace(" ", "").lower()

        # Retrieve variables from the NormDiff class, setting as attributes.
        if verbose:
            print("\tRetrieving Variables from NormDiff Instance...")
        try:
            # Ring radius and normalized power.
            self.rho_km_vals = NormDiff.rho_km_vals
            self.p_norm_vals = NormDiff.p_norm_vals
            
            # Phase of signal, negating do to mathematical conventions.
            self.phase_rad_vals = -NormDiff.phase_rad_vals

            # Ring opening angle.
            self.B_rad_vals = NormDiff.B_rad_vals

            # Spacecraft-to-Ring Intercept Point (RIP) distance.
            self.D_km_vals = NormDiff.D_km_vals

            # Ring azimuth angle.
            self.phi_rad_vals = NormDiff.phi_rad_vals

            # Frequency from the recieved signal.
            self.f_sky_hz_vals = NormDiff.f_sky_hz_vals

            # RIP velocity.
            self.rho_dot_kms_vals = NormDiff.rho_dot_kms_vals

            # Retrieve time variables (Earth, Ring, and Spacecraft ET).
            self.t_oet_spm_vals = NormDiff.t_oet_spm_vals  
            self.t_ret_spm_vals = NormDiff.t_ret_spm_vals         
            self.t_set_spm_vals = NormDiff.t_set_spm_vals  

            # Pole corrections in ring radius.       
            self.rho_corr_pole_km_vals = NormDiff.rho_corr_pole_km_vals

            # Timing corrections in ring radius.
            self.rho_corr_timing_km_vals = NormDiff.rho_corr_timing_km_vals

            # Ring longitude angle.
            self.phi_rl_rad_vals = NormDiff.phi_rl_rad_vals

            # Optical depth of diffraction profile.
            self.raw_tau_threshold_vals = NormDiff.raw_tau_threshold_vals

            # History from the NormDiff instance.
            self.dathist = NormDiff.history
        except AttributeError as errmes:
            raise AttributeError("NormDiff missing an attribute. %s" % errmes)

        # Compute mu: Sin(|B|)
        self.mu_vals = np.sin(np.abs(self.B_rad_vals))

        if verbose:
            print("\tCheck Variables for Errors...")
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

        if verbose:
            print("\tComputing Necessary Variables...")
        # Compute wavelength (km).
        self.lambda_sky_km_vals = SPEED_OF_LIGHT_KM / self.f_sky_hz_vals

        # Compute sampling distance (km)
        self.dx_km = self.rho_km_vals[1] - self.rho_km_vals[0]

        # Compute the complex transmittance.
        self.T_hat_vals = np.sqrt(np.abs(self.p_norm_vals))*np.exp(
            1j*self.phase_rad_vals)

        # Compute geometric qunatities for the Fresnel Scale.
        cb = np.cos(self.B_rad_vals)
        sb = np.sin(self.B_rad_vals)
        sp = np.sin(self.phi_rad_vals)

        # Compute the Fresnel Scale (km).
        self.F_km_vals = np.sqrt(0.5 * self.lambda_sky_km_vals *
                                 self.D_km_vals * (1 - cb*cb*sp*sp)/(sb*sb))
        
        # Compute the Normalized Equaivalent Width (See MTR86 Equation 20)
        self.norm_eq            = self.__func_dict[wtype]["normeq"]

        # Compute the window width. (See MTR86 Equations 19, 32, and 33).
        if bfac:
            omega = 2.0 * np.pi * self.f_sky_hz_vals
            alpha = (omega*omega) * (sigma*sigma)/(2.0 * self.rho_dot_kms_vals)
            P     = res / (alpha * (self.F_km_vals*self.F_km_vals))
            # The inverse exists only if P>1.
            if (np.min(P) <= 1.0):
                raise ValueError(
                    "\nWarning: Bad Points!\n\
                    \rEither rho_dot_kms_vals, F_km_vals, or res_km is\n\
                    \rtoo small, or sigma is too large. Request coarser\n\
                    \rresolution or set bfac=False as a keyword."
                )
            self.w_km_vals = self.norm_eq*np.abs(((P-1)*lambertw(
                np.exp(P/(1-P))*P/(1-P))+P)/(P-1)/alpha)
        else:
            self.w_km_vals = 2.0*self.norm_eq*self.F_km_vals*self.F_km_vals/res

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
        else:
            self.rng = np.array([np.min(rng),np.max(rng)])

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
        if self.verbose:
            print("\tRunning Fresnel Inversion...")
        self.T_vals = self.__ftrans(fwd=False)
        if self.verbose:
            sys.stdout.write("\x1b[2K")
            print("\tInversion Complete.")
        if self.fwd:
            if self.verbose:
                print("\tComputing Forward Transform...")
            self.T_hat_fwd_vals  = self.__ftrans(fwd=True)
            self.p_norm_fwd_vals = np.abs(
                self.T_hat_fwd_vals*self.T_hat_fwd_vals
            )
            self.phase_fwd_vals = -np.arctan2(np.imag(self.T_hat_fwd_vals),
                                              np.real(self.T_hat_fwd_vals))
            if self.verbose:
                sys.stdout.write("\x1b[2K")
                print("Forward Transform Complete.")
        if self.verbose:
            print("\tComputing Power and Phase...")

        # Compute power and phase.
        self.power_vals = np.abs(self.T_vals*self.T_vals)
        self.phase_vals = -np.arctan2(np.imag(self.T_vals),np.real(self.T_vals))

        # Compute regions of non-zero power.
        crange = (self.power_vals>0).nonzero()

        # Create empty array for normalized optical depth.
        tau = np.zeros(np.size(self.power_vals))

        # Compute the normalized optical depth.
        tau[crange] = -self.mu_vals[crange]*np.log(
            np.abs(self.power_vals[crange]))
        self.tau_vals = tau

        self.tau_threshold_vals = np.zeros(np.size(self.rho_km_vals))

        self.__trim_attributes(self.fwd)

        if self.verbose:
            print("\tWriting History...")

        input_vars = {
            "Diffraction Data": NormDiff,
            "Resolution (km)":  self.res,
        }

        input_kwds = {
            "Requested Range":      self.rngreq,
            "Requested Window":     self.wtype,
            "Use of Forward Model": self.fwd,
            "Use of Normalization": self.norm,
            "Use of Verbose":       self.verbose,
            "Use of b-factor":      self.bfac,
            "sigma value":          self.sigma,
            "Use of FFT":           self.fft,
            "Psi Approximation":    self.psitype
        }

        self.history = write_history_dict(input_vars, input_kwds, __file__)

        t2 = time.time()
        if self.verbose:
            print("\tDiffraction Correction Complete.")
        print("Computation Time: ",t2-t1,end="\r")

    def __rect(w_in, dx):
        """
            Rectangular Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute window function
        w_func = np.zeros(nw_pts) + 1.0
        return w_func

    def __coss(w_in, dx):
        """
            Squared Cosine Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Compute window function.
        w_func = np.cos(np.pi * x / w_in)**2
        return w_func

    def __kb20(w_in, dx):
        """
            Kaiser-Bessel 2.0 Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb20 is 2.0
        alpha  = 2.0 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb25(w_in, dx):
        """
            Kaiser-Bessel 2.5 Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb25 is 2.5.
        alpha  = 2.5 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb35(w_in, dx):
        """
            Kaiser-Bessel 3.5 Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kb35 is 3.5.
        alpha  = 3.5 * np.pi
        # Compute window function.
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kbmd20(w_in, dx):
        """
            Modified Kaiser-Bessel 2.0 Window Function.
        """
        # Window functions have an odd number of points.
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        # Compute argument of window function.
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        # Alpha value for kbmd20 is 2.0.
        alpha  = 2.0*np.pi
        # Compute window function.
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    def __kbmd25(w_in, dx):
        """
            Modified Kaiser-Bessel 2.5 window function.
        """
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
        "rect":   {"func": __rect,     "normeq": 1.00000000},
        "coss":   {"func": __coss,     "normeq": 1.50000000},
        "kb20":   {"func": __kb20,     "normeq": 1.49634231},
        "kb25":   {"func": __kb25,     "normeq": 1.65191895},
        "kb35":   {"func": __kb35,     "normeq": 1.92844639},
        "kbmd20": {"func": __kbmd20,   "normeq": 1.52048174},
        "kbmd25": {"func": __kbmd25,   "normeq": 1.65994218}
        }

    def __trim_attributes(self, fwd):
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

        # If the forward model was run, trim those attributes as well.
        if fwd:
            # Forward power
            self.p_norm_fwd_vals = self.p_norm_fwd_vals[crange]

            # Forward Transmittance and phase.
            self.T_hat_fwd_vals = self.T_hat_fwd_vals[crange]
            self.phase_fwd_vals = self.phase_fwd_vals[crange]

    def __fresinv(self, T_hat, ker, dx, f_scale):
        """
            Approximation Fresnel Inverse (MTR86 Equation 15)
        """
        T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
        return T

    def __fresinvfft(self, T_hat, ker, dx, f_scale):
        """
            FFT approximation of Fresnel Inverse (MTR86 Equation 15)
        """

        # Number of points in window.
        nw = np.size(T_hat)

        # FFT of data.
        fft_t_hat = np.fft.fft(T_hat)

        # FFT of weight Fresnel Kernel.
        fft_conv = np.fft.fft(ker)

        # Inverse FFT, assuming the convolution theorem is valid.
        inv_t_hat = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))

        # Scale factor outside of integral.
        inv_t_hat *= dx*(1.0+1.0j)/(2.0*f_scale)

        # Return midpoint value.
        T = inv_t_hat[int((nw-1)/2)+1]
        return T

    def __normalize(self, dx, ker, f_scale):
        """
            Window normalization function.
        """
        # Freespace Integral
        T1 = np.abs(np.sum(ker) * dx)

        # Normalization Factor
        norm_fact = np.sqrt(2.0) * f_scale / T1
        return norm_fact

    def __psi_func(self, kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw):
        # Compute Xi variable (MTR86 Equation 4b).
        xi = cbd * (r0 * cp0 - r * cp)

        # Compute Eta variable (MTR86 Equation 4c).
        eta = (r0*r0 + r*r - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

        psi_vals = kD * (np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi))
        return psi_vals

    def __psi_quartic(self, kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw):
        # Compute psi.
        psi_full = self.__psi_func(kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw)

        # Compute shifted ring radius and window width.
        x = (r-r0)
        w = np.max(x)-np.min(x)

        # Compute midpoints and endpoints.
        n1 = 0
        n2 = np.min((r0 - w/4 <= r).nonzero())
        n3 = np.max((r0 + w/4 >= r).nonzero())
        n4 = nw-1

        # Compute variables for psi polynomials.
        d_psi_half = psi_full[n3]-psi_full[n2]
        d_psi_full = psi_full[n4] - psi_full[n1]
        a_psi_half = (psi_full[n3]+psi_full[n2])/2
        a_psi_full = (psi_full[n1]+psi_full[n4])/2

        # Compute coefficients for the polynomial expansion.
        c1 = (8.0*d_psi_half-d_psi_full)/(3.0*w)
        c2 = 4.0*(16.0*a_psi_half-a_psi_full)/(3.0*w*w)
        c3 = 16.0*(d_psi_full-2.0*d_psi_half)/(3.0*w*w*w)
        c4 = 64.0*(a_psi_full-4.0*a_psi_half)/(3.0*w*w*w*w)

        # Compute quartic psi approximation.
        psi_vals = (c1 + c2*x)*x + (c3+c4*x)*x*x*x
        return psi_vals

    def __ftrans(self, fwd):
        # Retrieve Ring Radius.
        rho_km_vals = self.rho_km_vals

        # Retrieve ring azimuth angle.
        phi_rad_vals = self.phi_rad_vals

        # Retrieve window width, RIP distance, and Fresnel scale.
        w_km_vals = self.w_km_vals
        d_km_vals = self.D_km_vals
        F_km_vals = self.F_km_vals

        # Retrieve opening angle.
        B_rad_vals = self.B_rad_vals

        # Retrieve wavelength of transmitted signal
        lambda_km_vals = self.lambda_sky_km_vals

        # Retrieve FFT Boolean
        fft = self.fft

        # Retrieve psi approximation and verbose Boolean
        psitype = self.psitype
        verbose = self.verbose

        # Retrieve norm Boolean.
        norm = self.norm

        # Retrieve requested window type, starting point, and sample spacing.
        wtype = self.wtype
        start = self.start
        dx_km = self.dx_km

        # Retrieve number of points used.
        n_used = self.n_used

        # Compute product of wavenumber and RIP distance.
        kD_vals = 2.0 * np.pi * d_km_vals / lambda_km_vals

        # Compute Cosine of opening angle.
        cosb = np.cos(B_rad_vals)

        # Precompute variables for speed.
        cosbD = cosb/d_km_vals
        cosb2 = cosb*cosb

        # Compute cosine and sine of ring azimuth angle.
        cosphi0 = np.cos(phi_rad_vals)
        sinphi0 = np.sin(phi_rad_vals)

        # Precompute squares of variables for speed.
        dsq = d_km_vals*d_km_vals
        rsq = rho_km_vals*rho_km_vals

        # Define window function.
        fw = self.__func_dict[wtype]["func"]

        # Define normalization function and verbose message.
        nrm = self.__normalize
        mes = "Pt: %d  Tot: %d  Width: %d  Psi Iters: %d"

        # Set inverse function to FFT or Integration.
        if fft:
            finv = self.__fresinvfft
        else:
            finv = self.__fresinv
        # Set psi approximation.
        if (psitype == 'full'):
            psif = self.__psi_func
        elif (psitype == 'mtr4'):
            psif = self.__psi_quartic
        else:
            psif = self.__psi_func
        # If forward transform, adjust starting point by half a window.
        if fwd:
            T_in = self.T_vals
            T_out = T_in * 0.0
            w_max = np.max(w_km_vals[start:start + n_used])
            nw_fwd = int(np.ceil(w_max / (2.0 * dx_km)))
            start_fwd = int(start + nw_fwd)
            n_fwd = int(n_used - 2 * nw_fwd)
            start = start_fwd
            n_used = n_fwd
        else:
            T_in = self.T_hat_vals
            T_out = T_in * 0.0

        # Set the first computed point to the 'start' variable.
        center = start

        # Compute first window width and window function.
        w_init = w_km_vals[center]
        w_func = fw(w_init, dx_km)

        # Compute number of points in window function
        nw = np.size(w_func)

        # Computed range about the first point
        crange = np.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))-1

        # Compute first approximation for stationary phase.
        phi_s_rad = phi_rad_vals[center]

        # Compute Cosines and Sines of first approximation.
        cp = np.cos(phi_s_rad)
        sp = np.sin(phi_s_rad)

        # Factor used for first Newton-Raphson iteration
        dphi_fac = (cosb2*cosphi0*sinphi0/(1.0-cosb2*sinphi0*sinphi0))
        if (psitype == 'fresnel'):
            loop = 0
            for i in np.arange(n_used):
                # Current point being computed.
                center = start+i

                # Ring radius of current point.
                r0 = rho_km_vals[center]

                # Window width for current point.
                w = w_km_vals[center]

                # Window function for current point.
                w_func = fw(w, dx_km)

                # Number of points in current window.
                nw = np.size(w_func)

                # Computed range of points.
                crange = np.arange(int(center-(nw-1)/2),
                                   int(1+center+(nw-1)/2))-1

                # Computed ring radius range and Fresnel scale.
                r = rho_km_vals[crange]
                F = F_km_vals[center]

                # Compute psi for with stationary phase value
                psi_vals = (np.pi/2.0)*((r-r0)/F)*((r-r0)/F)

                # Compute kernel function for Fresnel inverse
                if fwd:
                    ker = w_func*np.exp(1j*psi_vals)
                else:
                    ker = w_func*np.exp(-1j*psi_vals)

                # Range of diffracted data that falls inside the window
                T = T_in[crange]

                # Compute 'approximate' Fresnel Inversion for current point
                T_out[center] = finv(T, ker, dx_km, F)

                # If normalization has been set, normalize the reconstruction
                if norm:
                    T_out[center] *= nrm(dx_km, ker, F)
                if verbose:
                    print(mes % (i, n_used, nw, loop), end="\r")
        else:
            for i in np.arange(n_used):
                # Current point being computed.
                center = start+i

                # Compute current radius and RIP distance.
                r0 = rho_km_vals[center]
                d2 = dsq[center]

                # Compute square of ring radius.
                r02 = rsq[center]

                # Precomputed variables for speed.
                cbd = cosbD[center]
                cb2 = cosb2[center]
                cp0 = cosphi0[center]
                sp0 = sinphi0[center]

                # Compute product of wavenumber and RIP Distance.
                kD = kD_vals[center]

                # Current window width and Fresnel scale.
                w = w_km_vals[center]
                F = F_km_vals[center]

                # If the window width has changed, recompute variables.
                if (np.abs(w_init - w) >= 2.0 * dx_km):
                    # Reset w_init and recompute window function.
                    w_init = w
                    w_func = fw(w, dx_km)

                    # Reset number of window points
                    nw = np.size(w_func)

                    # Computed range for current point
                    crange = np.arange(int(center-(nw-1)/2),
                                       int(1+center+(nw-1)/2))

                    # Ajdust ring radius by dx_km.
                    r = rho_km_vals[crange]

                    # Compute square of ring radius.
                    r2 = rsq[crange]

                    # First iteration of Newton-Raphson.
                    dphi_s_rad = dphi_fac[center] * (r - r0) / r0

                    # Compute new angle.
                    phi_s_rad = phi_rad_vals[center] - dphi_s_rad

                    # Compute Sines and Cosines of new angle.
                    cp = np.cos(phi_s_rad)
                    sp = np.sin(phi_s_rad)

                    # Compute r*cos(B)*D for efficiency
                    rcbd = r * cbd

                # If the window width has not changed, perform Newton-Raphson.
                else:
                    # Adjust computed range by dx_km.
                    crange += 1

                    # Ajdust ring radius by dx_km.
                    r = rho_km_vals[crange]

                    # Compute r*cos(B)*D for efficiency
                    rcbd = r * cbd

                    # Compute square of ring radius.
                    r2 = rsq[crange]

                    # Compute Xi variable (MTR86 Equation 4b).
                    xi = cbd * (r0 * cp0 - r * cp)

                    # Compute Eta variable (MTR86 Equation 4c).
                    eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                    # Compute Xi variable (MTR86 Equation 4b).
                    xi = cbd * (r0 * cp0 - r * cp)

                    # Compute Eta variable (MTR86 Equation 4c).
                    eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                    # Compute intermediate variables for partial derivatives.
                    v1 = rcbd * sp
                    v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                    v3 = np.sqrt(1.0 + 2.0*xi + eta)
                    v4 = rcbd * cp
                    v5 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                    # Compute variables used for second partial derivative.
                    dphia = (2.0*v4 + v5)/(2.0 * v3)
                    dphib = v4 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v3*v3*v3))

                    # Compute First and Second Partial Derivatives of psi
                    psi_d1 = (2.0*v1 + v2) / (2.0 * v3) - v1
                    psi_d2 = dphia - dphib

                    # Compute Newton-Raphson perturbation
                    dphi_s_rad = -psi_d1 / psi_d2
                    phi_s_rad += dphi_s_rad

                    # Compute Cosines and Sines new angle
                    cp = np.cos(phi_s_rad)
                    sp = np.sin(phi_s_rad)

                loop = 0
                # Perform Newton-Raphson on phi.
                while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
                    # Compute Xi variable (MTR86 Equation 4b).
                    xi = cbd * (r0 * cp0 - r * cp)

                    # Compute Eta variable (MTR86 Equation 4c).
                    eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                    # Compute intermediate variables for partial derivatives.
                    v1 = rcbd * sp
                    v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                    v3 = np.sqrt(1.0 + 2.0*xi + eta)
                    v4 = rcbd * cp
                    v5 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                    # Compute variables used for second partial derivative.
                    dphia = (2.0*v4 + v5)/(2.0 * v3)
                    dphib = v4 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v3*v3*v3))

                    # Compute First and Second Partial Derivatives of psi
                    psi_d1 = (2.0*v1 + v2) / (2.0 * v3) - v1
                    psi_d2 = dphia - dphib

                    # Compute Newton-Raphson perturbation
                    dphi_s_rad = -psi_d1 / psi_d2
                    phi_s_rad += dphi_s_rad

                    # Compute Cosines and Sines new angle
                    cp = np.cos(phi_s_rad)
                    sp = np.sin(phi_s_rad)

                    # Add one to loop variable for each iteration
                    loop += 1
                    if loop > 5:
                        break

                # Compute psi for with stationary phase value
                psi_vals = psif(kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw)

                # Compute kernel function for Fresnel inverse
                if fwd:
                    ker = w_func*np.exp(1j*psi_vals)
                else:
                    ker = w_func*np.exp(-1j*psi_vals)

                # Range of diffracted data that falls inside the window
                T = T_in[crange]

                # Compute 'approximate' Fresnel Inversion for current point
                T_out[center] = finv(T, ker, dx_km, F)

                # If normalization has been set, normalize the reconstruction
                if norm:
                    T_out[center] *= nrm(dx_km, ker, F)
                if verbose:
                    print(mes % (i, n_used, nw, loop), end="\r")
        return T_out
