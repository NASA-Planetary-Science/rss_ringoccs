# Import dependencies for the diffcorr module
import time
import numpy as np
from scipy.special import lambertw, iv
from scipy.constants import speed_of_light
from rss_ringoccs.tools.write_history_dict import write_history_dict

# Dictionary containing regions of interest within the Saturnian Rings.
region_dict = {
    'all':             [65000.0, 145000.0],
    'cringripples':    [77690.0, 77760.0],
    'encke':           [132900.0, 134200.0],
    'enckegap':        [132900.0, 134200.0],
    'janus':           [96200.0, 96800.0],
    'janusepimetheus': [96200.0, 96800.0],
    'maxwell':         [87410.0, 87610.0],
    'maxwellringlet':  [87410.0, 87610.0],
    'titan':           [77870.0, 77930.0],
    'titanringlet':    [77870.0, 77930.0],
    'huygens':         [117650.0, 117950.0],
    'huygensringlet':  [117650.0, 117950.0]
    }


class DiffractionCorrection(object):

    def __init__(self, NormDiff, res, rng="all", wtype="kb25", fwd=False,
                 norm=True, verbose=False, bfac=True, sigma=2.e-13,
                 fft=False, psitype="full"):
        t1 = time.time()

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

        # Set forward power variable to None in case it isn't defined later.
        self.p_norm_fwd_vals = None

        # Set forward phase and power to None as well.
        self.T_hat_fwd_vals = None
        self.phase_fwd_vals = None

        # Assign keywords/arguments to their corresponding attribute.
        self.res        = res
        self.wtype      = wtype.replace(" ", "").lower()
        self.rngreq     = rng
        self.norm       = norm
        self.verbose    = verbose
        self.fwd        = fwd
        self.fft        = fft
        self.bfac       = bfac
        self.sigma      = sigma
        self.psitype    = psitype.replace(" ", "").lower()

        # Retrieve variables from the NormDiff class, setting as attributes.
        try:
            self.rho_km_vals                = NormDiff.rho_km_vals
            self.p_norm_vals                = NormDiff.p_norm_vals
            self.phase_rad_vals             = -NormDiff.phase_rad_vals
            self.B_rad_vals		            = NormDiff.B_rad_vals
            self.D_km_vals		            = NormDiff.D_km_vals
            self.f_sky_hz_vals              = NormDiff.f_sky_hz_vals
            self.phi_rad_vals               = NormDiff.phi_rad_vals
            self.rho_dot_kms_vals           = NormDiff.rho_dot_kms_vals
            self.t_oet_spm_vals             = NormDiff.t_oet_spm_vals  
            self.t_ret_spm_vals             = NormDiff.t_ret_spm_vals         
            self.t_set_spm_vals             = NormDiff.t_set_spm_vals         
            self.rho_corr_pole_km_vals      = NormDiff.rho_corr_pole_km_vals
            self.rho_corr_timing_km_vals    = NormDiff.rho_corr_timing_km_vals
            self.phi_rl_rad_vals            = NormDiff.phi_rl_rad_vals
            self.raw_tau_threshold_vals     = NormDiff.raw_tau_threshold_vals
            self.dathist                    = NormDiff.history
        except AttributeError as errmes:
            raise AttributeError(
                "NormDiff missing an attribute. %s" % errmes
            )

        # Compute mu: Sin(|B|)
        self.mu_vals = np.sin(np.abs(self.B_rad_vals))

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
        self.T_hat_vals         = np.sqrt(np.abs(self.p_norm_vals))*np.exp(
            1j*self.phase_rad_vals)

        # Compute the Fresnel Scale (See MTR86 Equation 6)
        cb                      = np.cos(self.B_rad_vals)
        sb                      = np.sin(self.B_rad_vals)
        sp                      = np.sin(self.phi_rad_vals)
        self.F_km_vals          = np.sqrt(0.5 * self.lambda_sky_km_vals *
            self.D_km_vals * (1 - (cb*cb) * (sp*sp)) / (sb*sb))
        
        # Compute the Normalized Equaivalent Width (See MTR86 Equation 20)
        self.norm_eq            = self.__func_dict[wtype]["normeq"]

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
        self.T_vals = self.__ftrans(fwd=False)
        if self.fwd:
            self.T_hat_fwd_vals  = self.__ftrans(fwd=True)
            self.p_norm_fwd_vals = np.abs(
                self.T_hat_fwd_vals*self.T_hat_fwd_vals
            )
            self.phase_fwd_vals = -np.arctan2(np.imag(self.T_hat_fwd_vals),
                                              np.real(self.T_hat_fwd_vals))

        self.power_vals = np.abs(self.T_vals*self.T_vals)
        self.phase_vals = -np.arctan2(np.imag(self.T_vals),np.real(self.T_vals))
        crange          = (self.power_vals>0).nonzero()
        tau             = np.zeros(np.size(self.power_vals))
        tau[crange]     = -self.mu_vals[crange]*np.log(
            np.abs(self.power_vals[crange]))
        self.tau_vals   = tau

        self.__trim_attributes(self.fwd)

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
        "rect"    : {"func" : __rect,     "normeq" : 1.00000000},
        "coss"    : {"func" : __coss,     "normeq" : 1.50000000},
        "kb20"    : {"func" : __kb20,     "normeq" : 1.49634231},
        "kb25"    : {"func" : __kb25,     "normeq" : 1.65191895},
        "kb35"    : {"func" : __kb35,     "normeq" : 1.92844639},
        "kbmd20"  : {"func" : __kbmd20,   "normeq" : 1.52048174},
        "kbmd25"  : {"func" : __kbmd25,   "normeq" : 1.65994218}
        }

    def __trim_attributes(self,fwd):
        # Get rid of uncomputed values and keep only what was processed.
        start  = self.start                 # Starting Point.
        n_used = self.n_used                # Number of point used.
        crange = np.arange(n_used)+start    # Processed range.

        self.rho_km_vals                = self.rho_km_vals[crange]
        self.p_norm_vals                = self.p_norm_vals[crange]
        self.phase_rad_vals             = self.phase_rad_vals[crange]
        self.B_rad_vals                 = self.B_rad_vals[crange]
        self.D_km_vals                  = self.D_km_vals[crange]
        self.f_sky_hz_vals              = self.f_sky_hz_vals[crange]
        self.phi_rad_vals               = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals           = self.rho_dot_kms_vals[crange]
        self.T_hat_vals                 = self.T_hat_vals[crange]
        self.F_km_vals                  = self.F_km_vals[crange]
        self.w_km_vals                  = self.w_km_vals[crange]
        self.mu_vals                    = self.mu_vals[crange]
        self.lambda_sky_km_vals         = self.lambda_sky_km_vals[crange]
        self.T_vals                     = self.T_vals[crange]
        self.power_vals                 = self.power_vals[crange]
        self.tau_vals                   = self.tau_vals[crange]
        self.phase_vals                 = self.phase_vals[crange]
        self.t_oet_spm_vals             = self.t_oet_spm_vals[crange]
        self.t_ret_spm_vals             = self.t_ret_spm_vals[crange]
        self.t_set_spm_vals             = self.t_set_spm_vals[crange]
        self.rho_corr_pole_km_vals      = self.rho_corr_pole_km_vals[crange]
        self.rho_corr_timing_km_vals    = self.rho_corr_timing_km_vals[crange]
        self.phi_rl_rad_vals            = self.phi_rl_rad_vals[crange]
        self.raw_tau_threshold_vals     = self.raw_tau_threshold_vals[crange]

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
        nw              = np.size(T_hat)
        fft_t_hat       = np.fft.fft(T_hat)
        fft_conv        = np.fft.fft(ker)
        inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
        inv_t_hat      *= dx*(1.0+1.0j)/(2.0*f_scale)
        T               = inv_t_hat[int((nw-1)/2)+1]
        return T

    # Window Normalization Function.
    def __normalize(self,dx,ker,f_scale):
        T1        = np.abs(np.sum(ker) * dx)   # Freespace Integral
        norm_fact = np.sqrt(2.0) * f_scale / T1         # Normalization Factor
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

    def __psi_func_fwd(self, kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw):
        # Compute Xi variable (MTR86 Equation 4b).
        xi = cbd * (r0 * cp0 - r * cp)

        # Compute Eta variable (MTR86 Equation 4c).
        eta = (r0*r0 + r*r - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

        psi_vals = -kD * (np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi))
        return psi_vals

    def __psi_quartic_fwd(self, kD, r, r0, cbd, d2, cp0, sp0, cp, sp, w, nw):
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
        psi_vals *= -1.0
        return psi_vals

    def __ftrans(self,fwd):
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

        # Define normalization function.
        nrm = self.__normalize

        # Set inverse function to FFT or Integration.
        if fft:
            finv = self.__fresinvfft
        else:
            finv = self.__fresinv

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
            # Set psi approximation.
            if (psitype == 'full'):
                psif = self.__psi_func_fwd
            elif (psitype == 'mtr4'):
                psif = self.__psi_quartic_fwd
            else:
                psif = self.__psi_func_fwd
            mes = "Pt: %d  Tot: %d  Width: %d  Psi Iters: %d Forward     "
        else:
            T_in = self.T_hat_vals
            T_out = T_in * 0.0
            # Set psi approximation.
            if (psitype == 'full'):
                psif = self.__psi_func
            elif (psitype == 'mtr4'):
                psif = self.__psi_quartic
            else:
                psif = self.__psi_func
            mes = "Pt: %d  Tot: %d  Width: %d  Psi Iters: %d Inversion   "

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
                ker = w_func*np.exp(-1j*psi_vals)

                # Range of diffracted data that falls inside the window
                T = T_in[crange]

                # Compute 'approximate' Fresnel Inversion for current point
                T_out[center] = finv(T_hat, ker, dx_km, F)

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

                # If the window width has not changed, perform Newton-Raphson.
                else:
                    # Adjust computed range by dx_km.
                    crange += 1

                    # Ajdust ring radius by dx_km.
                    r = rho_km_vals[crange]

                    # Compute square of ring radius.
                    r2 = rsq[crange]

                    # Compute Xi variable (MTR86 Equation 4b).
                    xi = cbd * (r0 * cp0 - r * cp)

                    # Compute Eta variable (MTR86 Equation 4c).
                    eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                    # Compute intermediate variables for partial derivatives.
                    v1 = r * cbd * sp
                    v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                    v4 = np.sqrt(1.0 + 2.0*xi + eta)
                    v5 = cbd * r * cp
                    v6 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                    # Compute variables used for second partial derivative.
                    dphia = (2.0*v5 + v6)/(2.0 * v4)
                    dphib = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))

                    # Compute First and Second Partial Derivatives of psi
                    psi_d1 = (2.0*v1 + v2) / (2.0 * v4) - v1
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
                    v1 = r * cbd * sp
                    v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                    v4 = np.sqrt(1.0 + 2.0*xi + eta)
                    v5 = cbd * r * cp
                    v6 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                    # Compute variables used for second partial derivative.
                    dphia = (2.0*v5 + v6)/(2.0 * v4)
                    dphib = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))

                    # Compute First and Second Partial Derivatives of psi
                    psi_d1 = (2.0*v1 + v2) / (2.0 * v4) - v1
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
