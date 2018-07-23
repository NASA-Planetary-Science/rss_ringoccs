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
    
    Window (Taper) Functions:
        rect................Rectangular window.
        coss................Squared cossine window.
        kb20................Kaiser-Bessel 2.0 window.
        kb25................Kaiser-Bessel 2.5 window.
        kb35................Kaiser-Bessel 3.5 window.
        kbmd20..............Modified Kaiser-Bessel 2.0 window.
        kbmd25..............Modified Kaiser-Bessel 2.5 window.
    
    Special Functions:
        fresnel_sin.........The Fresnel sine integral.
        fresnel_cos.........The Fresnel cosine integral.
        sq_well_solve.......Diffraction pattern through square well.

    Mathematical Functions:
        compute_norm_eq.....Computes the normalized equivalent width.
        get_norm_eq.........Quickly retrieve pre-computed normalized
                            equivalent widths from strings with the
                            name of common window functions.
        resolution_inverse..Computes the inverse of the function
                            y = x/(exp(-x)+x-1)
        power_func..........Compute power from complex transmittance.
        phase_func..........Compute phase from complex transmittance.
        tau_func............Compute normalized optical depth from the
                            complex transmittance.
        wker................Computes a weighted kernel function.
        freq_wav............Convert frequency to wavelength, and
                            vice-versa. Kilometers or Herts only.
        fresnel_scale.......Compute the Fresnel scale.
    
    Miscellaneous Functions:
        get_range_request...Computes an array of the form [a,b] from
                            a given array, list, or from a set of
                            allowed strings.
        get_range_actual....Given an array of numbers (usually the
                            radial range of the data), a range
                            request, and a window width, compute the
                            allowed range of processing.
"""
# Import dependencies for the diffcorr module
import time,os,sys,platform
import numpy as np
from scipy.special import lambertw, iv
from scipy import interpolate
from scipy.constants import speed_of_light
sys.path.append("%s/../../" % os.path.dirname(os.path.realpath(__file__)))
from rss_ringoccs.tools import check_boole,check_real
from rss_ringoccs.tools import check_pos_real,check_complex
import pdb

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
        Class: rec_data
        Purpose: This object contains all of input and output variables from
        diffraction reconstruction, including geometry, data, and reconstructed
        data.
        Attributes:
            rho_km_vals (np.ndarray): Ring intercept radius values, in km,
                at final spacing specified in dr_km
            p_norm_vals (np.ndarray): Power normalized to 1. This is the diffraction
                pattern that is input to the Fresnel inversion step
            phase_rad_vals (np.ndarray): Phase of the complex signal, in radians.
                This is the other part of the diffraction pattern that is input to
                the Fresnel Inversion step
            B_rad_vals (np.ndarray): Ring opening angle of Saturn
            D_km_vals (np.ndarray): Spacecraft-RIP (Ring Intercept Point) distance
            f_sky_hz_vals (np.ndarray): Predicted sky frequency, in Hz.
            phi_rad_vals (np.ndarray): Observed ring azimuth, in radians
            rho_dot_kms_vals (np.ndarray): Ring intercept point velocity
	"""

    def __init__(self, NormDiff,res,wtype,bfac=True,sigma=False):
        """
            Arguments:
                NormDiff: Instance of the class NormDiff containing 
        """

        if not check_pos_real(res):
            raise TypeError("res must be a positive number")
        if not (type(wtype)==type("Hi!")):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        if not check_boole(bfac):
            raise TypeError("bfac must be Boolean: True/False")
        if (sigma == False): pass
        else:
            if (not check_pos_real(sigma)):
                raise TypeError("sigma must be a positive number.")
            else: pass

        self.res                     = None
        self.wtype                   = None
        self.rho_km_vals             = None
        self.p_norm_vals             = None
        self.phase_rad_vals          = None
        self.B_rad_vals              = None
        self.D_km_vals               = None
        self.f_sky_hz_vals           = None
        self.phi_rad_vals            = None
        self.rho_dot_kms_vals        = None
        self.T_hat_vals              = None
        self.F_km_vals               = None
        self.w_km_vals               = None
        self.mu_vals                 = None
        self.lambda_sky_km_vals      = None
        self.dx_km                   = None
        self.norm_eq                 = None
        self.t_oet_spm_vals          = None
        self.t_ret_spm_vals          = None
        self.t_set_spm_vals          = None
        self.rho_corr_pole_km_vals   = None
        self.rho_corr_timing_km_vals = None
        self.phi_rl_rad_vals         = None

        self.res   = res
        self.wtype = wtype.replace(" ", "").lower()

        self.rho_km_vals      = NormDiff.rho_km_vals
        self.p_norm_vals      = NormDiff.p_norm_vals
        self.phase_rad_vals   = -NormDiff.phase_rad_vals
        self.B_rad_vals		  = NormDiff.B_rad_vals
        self.D_km_vals		  = NormDiff.D_km_vals
        self.f_sky_hz_vals    = NormDiff.f_sky_hz_vals
        self.phi_rad_vals     = NormDiff.phi_rad_vals
        self.rho_dot_kms_vals = NormDiff.rho_dot_kms_vals
        self.mu_vals          = np.sin(np.abs(self.B_rad_vals))
        try:
            self.t_oet_spm_vals          = NormDiff.t_oet_spm_vals         
            self.t_ret_spm_vals          = NormDiff.t_ret_spm_vals         
            self.t_set_spm_vals          = NormDiff.t_set_spm_vals         
            self.rho_corr_pole_km_vals   = NormDiff.rho_corr_pole_km_vals
            self.rho_corr_timing_km_vals = NormDiff.rho_corr_timing_km_vals
            self.phi_rl_rad_vals         = NormDiff.phi_rl_rad_vals
        except AttributeError:
            try:
                self.t_oet_spm_vals          = NormDiff.spm_vals
                self.t_ret_spm_vals          = NormDiff.t_ret_spm_vals         
                self.t_set_spm_vals          = NormDiff.t_set_spm_vals         
                self.rho_corr_pole_km_vals   = NormDiff.rho_corr_pole_km_vals
                self.rho_corr_timing_km_vals = NormDiff.rho_corr_timing_km_vals
                self.phi_rl_rad_vals         = NormDiff.phi_rl_rad_vals
            except AttributeError:
                pass

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

        self.lambda_sky_km_vals = speed_of_light*0.001 / self.f_sky_hz_vals
        self.dx_km              = self.rho_km_vals[1] - self.rho_km_vals[0]
        self.T_hat_vals         = np.sqrt(np.abs(self.p_norm_vals))*np.exp(1j*self.phase_rad_vals)
        cb                      = np.cos(self.B_rad_vals)
        sb                      = np.sin(self.B_rad_vals)
        sp                      = np.sin(self.phi_rad_vals)
        self.F_km_vals          = np.sqrt(0.5 * self.lambda_sky_km_vals *\
            self.D_km_vals * (1 - (cb*cb) * (sp*sp)) / (sb*sb))
        self.norm_eq            = self.__func_dict[wtype]

        if bfac:
            if (not sigma):
                sigma = 2.e-13
            omega = 2.0 * np.pi * self.f_sky_hz_vals
            alpha = (omega*omega) * (sigma*sigma) / (2.0 * self.rho_dot_kms_vals)
            P     = res / (alpha * (self.F_km_vals*self.F_km_vals))
            # The inverse exists only if P>1.
            if (np.min(P) < 1.0001):
                raise ValueError("Warning: Bad Points! Either rho_dot_kms_vals,\n\
                F_km_vals, or res_km is to small. Exclude\n\
                points or request coarser resolution.")
            self.w_km_vals = self.norm_eq*np.abs(((P-1)*lambertw(np.exp(P/(1-P))*P/(1-P))+P)/(P-1)/alpha)
        else:
            self.w_km_vals = 2.0*self.norm_eq*self.F_km_vals*self.F_km_vals/res

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

        if (type(rng) == type('Hello')):
            reg = rng.replace(" ","").lower()
            if (reg in region_dict):
                self.rng = np.array(region_dict[reg])
            else:
                print("Illegal range request.\nAllowed strings are:")
                for key in region_dict: print(key)
                raise ValueError("Illegal String. See list for valid strings.")
        elif (np.size(self.rngreq) < 2):
            raise TypeError("Only one value given for range. Use rng = [a,b]")
        else:
            self.rng = np.array([np.min(self.rngreq),np.max(self.rngreq)])

        # Compute the starting point and the number of points used.
        rho         = self.rho_km_vals          # Ring radius.
        w_max       = np.max(self.w_km_vals)    # Largest window used.

        # Compute the smallest and largest allowed radii for reconstruction.
        rho_min_lim = np.min(rho)+np.ceil(w_max/2.0)
        rho_max_lim = np.max(rho)-np.ceil(w_max/2.0)

        # Compute the smallest and largest values within requested range.
        rho_start   = rho[np.min((rho >= np.min(self.rng)).nonzero())]
        rho_end     = rho[np.max((rho <= np.max(self.rng)).nonzero())]

        # Compute the start and endpoint of reconstruction.
        rho_min     = np.max([rho_min_lim,rho_start])
        rho_max     = np.min([rho_max_lim,rho_end])
        self.start  = int(np.min((rho >= rho_min).nonzero()))
        self.finish = int(np.max((rho <= rho_max).nonzero()))
        self.n_used = 1 + (self.finish - self.start)

        # Delete unnecessary variables for clarity and memory.
        del rho,w_max,rho_min_lim,rho_max_lim,rho_start,rho_end,rho_min,rho_max

        #self.__trim_inputs()

        if verbose: self.T_vals = self.__finv_v()
        else: self.T_vals = self.__finv()
        if fwd:
            if norm: self.T_hat_fwd_vals = self.__ffwd_n_v()
            else:    self.T_hat_fwd_vals = self.__ffwd_v()
            self.p_norm_fwd_vals = power_func(self.T_hat_fwd_vals)
            self.phase_fwd_vals  = phase_func(self.T_hat_fwd_vals)

        self.power_vals = np.abs(self.T_vals*self.T_vals)
        self.phase_vals = -np.arctan2(np.imag(self.T_vals),np.real(self.T_vals))
        crange          = (self.power_vals>0).nonzero()
        tau             = np.zeros(np.size(self.power_vals))
        tau[crange]     = -self.mu_vals[crange] * np.log(np.abs(self.power_vals[crange]))
        self.tau_vals   = tau

        self.__trim_attributes(fwd)

        del fwd, verbose, psitype, norm

        self.history = self.__write_hist_dict()

        t2 = time.time()
        sys.stdout.write("\033[K")
        print("Computation Time: ",t2-t1,end="\r")

    def __rect(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        w_func = np.zeros(nw_pts) + 1.0
        return w_func

    def __coss(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        w_func = np.cos(np.pi * x / w_in)**2
        return w_func

    def __kb20(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.0*np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb25(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5*np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb35(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 3.5 * np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kbmd20(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.0*np.pi
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    def __kbmd25(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5 * np.pi
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    def __trim_attributes(self,fwd):
        start  = self.start
        n_used = self.n_used
        crange = np.arange(n_used)+start

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

    def __fresinv(self,T_hat,ker,dx,f_scale):
        T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
        return T

    def __fresinvfft(self,T_hat,ker,dx,f_scale):
        nw = np.size(T_hat)
        fft_t_hat       = np.fft.fft(T_hat)
        fft_conv        = np.fft.fft(ker)
        inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
        inv_t_hat      *= dx*(np.complex(1.0,1.0))/(2.0*f_scale)
        T               = inv_t_hat[int((nw-1)/2)]
        return T

    def __normalize(self,r,w_func,f_scale):
        x         = r-np.mean(r)
        drho      = r[1]-r[0]
        f_scale   = f_scale
        psi       = (np.pi / 2.0) * ((x / f_scale)*(x / f_scale))
        ker       = np.exp(-1j * psi)
        T1        = np.abs(np.sum(w_func * ker) * drho)
        norm_fact = np.sqrt(2.0) * f_scale / T1
        return norm_fact
    
    def __normalize_do_nothing(self,x,y,z):
        return 1

    __func_dict = {
        "__rect" : {"func" : __rect, "normeq" : 1.00000000},
        "__coss" : {"func" : __coss, "normeq" : 1.50000000},
        "__kb20" : {"func" : __kb20, "normeq" : 1.49634231},
        "__kb25" : {"func" : __kb25, "normeq" : 1.65191895},
        "__kb35" : {"func" : __kb35, "normeq" : 1.92844639},
        "__kbmd20" : {"func" : __kbmd20, "normeq" : 1.52048174},
        "__kbmd25" : {"func" : __kbmd25, "normeq" : 1.65994218}
        }

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
        # Calculate the corrected complex amplitude, point by point
        T_vals     = T_hat_vals * 0.0
        w_init     = w_vals[start]
        w_func     = fw(w_init,dx)
        nw         = np.size(w_func)
        crange     = np.arange(int(start-(nw-1)/2),int(1+start+(nw-1)/2))-1
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d2     = dsq[center]
            cb_d   = cosb_D[center]
            cb2    = cosb2[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (w_init - w>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = ((cb2)*cp0*sp0/(1.0-(cb2) * (sp0*sp0))) * (r - r0) / r0
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange    +=1
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb_d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb_d * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb_d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb_d * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            phi_s_rad1 = phi_s_rad
            
            # Compute psi and then compute the forward model.
            xi       = (cb_d) * (r0*cp0 - r*cp)
            eta      = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / d2
            psi_vals = kD * (np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi))
            F        = F_vals[center]
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center]   = finv(T_hat,ker,dx,F)
            T_vals[center]  *= nrm(r,w_func,F)
        return T_vals

    def __finv_v(self):
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
        # Calculate the corrected complex amplitude, point by point
        T_vals     = T_hat_vals * 0.0
        w_init     = w_vals[start]
        w_func     = fw(w_init,dx)
        nw         = np.size(w_func)
        crange     = np.arange(int(start-(nw-1)/2),int(1+start+(nw-1)/2))-1
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d2     = dsq[center]
            cb_d   = cosb_D[center]
            cb2    = cosb2[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (w_init - w>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = ((cb2)*cp0*sp0/(1.0-(cb2) * (sp0*sp0))) * (r - r0) / r0
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange    +=1
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb_d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb_d * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb_d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb_d * sp
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb_d * cp
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            phi_s_rad1 = phi_s_rad
            
            # Compute psi and then compute the forward model.
            xi       = (cb_d) * (r0*cp0 - r*cp)
            eta      = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / d2
            psi_vals = kD * (np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi))
            F        = F_vals[center]
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center]   = finv(T_hat,ker,dx,F)
            T_vals[center]  *= nrm(r,w_func,F)
            print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
            % (i,n_used,nw,loop),end="\r")
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