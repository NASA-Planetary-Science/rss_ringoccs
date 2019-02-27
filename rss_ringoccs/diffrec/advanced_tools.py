import numpy as np
from scipy import interpolate
from . import diffraction_correction, special_functions, window_functions
from rss_ringoccs.tools import CSV_tools, error_check


class CompareTau(object):
    def __init__(self, geo, cal, dlp, tau, res, rng='all', wtype="kbmd20",
                 fwd=False, bfac=True, sigma=2.e-13, verbose=False, norm=True,
                 psitype='fresnel4', res_factor=0.75):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.CompareTau"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(cal, str, "cal", fname)
        error_check.check_type(dlp, str, "dlp", fname)
        error_check.check_type(tau, str, "tau", fname)
        error_check.check_type(fwd, bool, "fwd", fname)
        error_check.check_type(norm, bool, "norm", fname)
        error_check.check_type(bfac, bool, "bfac", fname)
        error_check.check_type(wtype, str, "wtype", fname)
        error_check.check_type(psitype, str, "psitype", fname)

        # Try converting certain inputs into floating point numbers.
        res = error_check.check_type_and_convert(res, float, "res", fname)
        sigma = error_check.check_type_and_convert(sigma, float, "sigma", fname)
        res_factor = error_check.check_type_and_convert(res_factor, float,
                                                        "res_factor", fname)

        # Check that these variables are positive.
        error_check.check_positive(res, "res", fname)
        error_check.check_positive(sigma, "sigma", fname)
        error_check.check_positive(res_factor, "res_factor", fname)

        # Check that the requested window type is a legal input.
        if not isinstance(wtype, str):
            erm = ""
            for key in func_dict:
                erm = "%s\t\t'%s'\n" % (erm, key)
            raise TypeError(
                "\n\tError Encountered:\n"
                "\t\trss_ringoccs.diffrec.CompareTau\n\n"
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
                    "\n\tError Encountered:\n"
                    "\t\trss_ringoccs.diffrec.CompareTau\n\n"
                    "\tIllegal string used for wtype.\n"
                    "\tYour string: '%s'\n"
                    "\tAllowed Strings:\n%s" % (wtype, erm)
                )
            else:
                pass

        # Check that range and psitype are legal inputs.
        rng = error_check.check_range_input(rng, fname)
        psitype = error_check.check_psitype(psitype, fname)

        data = CSV_tools.ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
        rec = diffraction_correction.DiffractionCorrection(
            data, res, rng=rng, wtype=wtype,
            fwd=fwd, norm=norm, bfac=bfac,
            sigma=sigma, psitype=psitype,
            res_factor=res_factor, verbose=verbose
        )

        self.rho_km_vals = rec.rho_km_vals
        tr = data.tau_rho
        rmin = np.min(self.rho_km_vals)
        rmax = np.max(self.rho_km_vals)
        tau_rstart = int(np.min((tr-rmin>=0).nonzero()))
        tau_rfin = int(np.max((rmax-tr>=0).nonzero()))

        self.rho_km_vals = rec.rho_km_vals
        self.power_vals = rec.power_vals
        self.phase_vals = rec.phase_vals
        self.tau_vals = rec.tau_vals
        self.tau_power = data.power_vals[tau_rstart:tau_rfin+1]
        self.tau_phase = data.phase_vals[tau_rstart:tau_rfin+1]
        self.tau_tau = data.tau_vals[tau_rstart:tau_rfin+1]
        self.wtype = wtype
        self.res = res


class FindOptimalResolution(object):
    def __init__(self, geo, cal, dlp, tau, sres, dres, nres, fwd=False,
                 norm=True, bfac=True, sigma=2.e-13, psitype="fresnel4",
                 rng='all', wlst=['kbmd20'], res_factor=0.75, verbose=True):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.FindOptimalResolution"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(cal, str, "cal", fname)
        error_check.check_type(dlp, str, "dlp", fname)
        error_check.check_type(tau, str, "tau", fname)
        error_check.check_type(fwd, bool, "fwd", fname)
        error_check.check_type(norm, bool, "norm", fname)
        error_check.check_type(bfac, bool, "bfac", fname)
        error_check.check_type(wlst, list, "wlst", fname)
        error_check.check_type(psitype, str, "psitype", fname)
        sres = error_check.check_type_and_convert(sres, float, "sres", fname)
        dres = error_check.check_type_and_convert(dres, float, "dres", fname)
        nres = error_check.check_type_and_convert(nres, int, "nres", fname)
        sigma = error_check.check_type_and_convert(sigma, float, "sigma", fname)
        res_factor = error_check.check_type_and_convert(res_factor, float, 
                                                        "res_factor", fname)

        error_check.check_positive(sres, "sres", fname)
        error_check.check_positive(dres, "dres", fname)
        error_check.check_positive(nres, "nres", fname)
        error_check.check_positive(sigma, "sigma", fname)

        # Check that the requested range is a legal input.
        rng = check_range_input(rng, fname)

        if (not all(isinstance(x, str) for x in wlst)):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\t\trss_ringoccs.diffrec.FindOptimalResolution\n\n"
                "\twlst must be a list of strings.\n"
            )
        else:
            for i in range(len(wlst)):
                wlst[i] = wlst[i].replace(" ", "").replace("'", "")
                wlst[i] = wlst[i].replace('"', "").lower()

            if not all((wtype in func_dict) for wtype in wlst):
                erm = ""
                for key in func_dict:
                    erm = "%s\t\t'%s'\n" % (erm, key)
                raise ValueError(
                    "\n\tError Encountered:\n"
                    "\t\trss_ringoccs.diffrec.FindOptimalResolution\n\n"
                    "\tIllegal string used for wtype.\n"
                    "\tYour strings: '%s'\n"
                    "\tAllowed Strings:\n%s" % (wtype, erm)
                )
            else:
                pass

        nwins = np.size(wlst)
        self.linf = np.zeros((nres, nwins))
        self.l2 = np.zeros((nres, nwins))
        self.ideal_res = np.zeros((nwins))
        eres = sres + (nres-1)*dres
        res = sres
        data = ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
        rec = diffraction_correction(data, res, rng=rng,
                                    wtype="kb35", verbose=False)
        start = int(np.min((data.tau_rho-np.min(rec.rho_km_vals)>=0).nonzero()))
        fin = int(np.max((np.max(rec.rho_km_vals)-data.tau_rho>=0).nonzero()))
        tau_power = data.power_vals[start:fin+1]
        for i in np.arange(nres):
            for j in range(nwins):
                wtype = wlst[j]
                recint = diffraction_correction(data, res, rng=rng, wtype=wtype)
                p_int = np.abs(recint.power_vals - tau_power)
                self.linf[i,j] = np.max(p_int)
                self.l2[i,j] = np.sqrt(np.sum(p_int*p_int)*recint.dx_km)

                if verbose:
                    print("%s %f %s %f %s %s"
                          % ('Res:',res,'Max:',eres,"WTYPE:",wtype))

            res += dres
        self.ideal_res = sres+dres*np.argmin(self.linf, axis=0)


class SquareWellFromGEO(object):
    def __init__(self, geo, lambda_km, res, rho, width, dx_km_desired=0.25,
                 occ="other", wtype='kb25', fwd=False, norm=True, bfac=True,
                 verbose=True, psitype='full'):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.SquareWellFromGEO"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(norm, bool, "norm", fname)
        error_check.check_type(bfac, bool, "bfac", fname)
        error_check.check_type(occ, str, "occ", fname)
        error_check.check_type(fwd, bool, "fwd", fname)
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(wtype, str, "wtype", fname)
        error_check.check_type(psitype, str, "psitype", fname)
        res = error_check.check_type_and_convert(res, float, "res", fname)
        width = error_check.check_type_and_convert(width, float, "width", fname)
        lambda_km = error_check.check_type_and_convert(lambda_km, float,
                                                       "lambda_km", fname)
        dx_km_desired = error_check.check_type_and_convert(dx_km_desired, float,
                                                           "dx_km_desired", fname)

        error_check.check_positive(res, "res", fname)
        error_check.check_positive(width, "width", fname)
        error_check.check_positive(dx_km_desired, "dx_km_desired", fname)
        error_check.check_positive(lambda_km, "lambda_km", fname)

        data = CSV_tools.get_geo(geo, verbose=verbose)
        occ = occ.replace(" ", "").replace("'", "").replace('"', "")
        occ = occ.lower()

        if verbose:
            print("Retrieving Variables...")

        # Retrieve rho and phi from geo data.
        self.rho_km_vals = np.array(data.rho_km_vals)
        self.phi_rad_vals = np.deg2rad(np.array(data.phi_ora_deg_vals))

        # Retrieve B and D from geo data.
        self.D_km_vals = np.array(data.D_km_vals)
        self.B_rad_vals = np.deg2rad(np.array(data.B_deg_vals))

        # Retrieve drho/dt.
        self.rho_dot_kms_vals = np.array(data.rho_dot_kms_vals)

        # Run error check on variables from GEO file.
        error_check.check_is_real(self.D_km_vals, "D_km_vals", fname)
        error_check.check_is_real(self.B_rad_vals, "B_rad_vals", fname)
        error_check.check_is_real(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_is_real(self.phi_rad_vals, "phi_rad_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals,
                                  "rho_dot_km_vals", fname)

        error_check.check_positive(self.D_km_vals, "D_km_vals", fname)
        error_check.check_positive(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_two_pi(self.B_rad_vals, "B_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rad_vals, "phi_rad_vals",
                                 fname, deg=False)

        if verbose:
            print("Computing Variables...")

        if (occ == 'ingress'):
            crange = (self.rho_dot_kms_vals < 0.0).nonzero()
        elif (occ == 'egress'):
            crange = (self.rho_dot_kms_vals > 0.0).nonzero()
        else:
            crange_e = (self.rho_dot_kms_vals > 0.0).nonzero()
            crange_i = (self.rho_dot_kms_vals < 0.0).nonzero()
            n_e      = np.size(crange_e)
            n_i      = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals has positive and negative values.
                        \r\tThis is likely a chord occultation.
                        \r\tSet occ='ingress' or occ='egress'.
                    """ % fname
                )
            elif (n_e == 0) and (n_i == 0):
                raise ValueError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\trho_dot_kms_vals is either zero or empty.
                    """ % fname
                )
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else:
                raise TypeError(
                    """
                        \r\tError Encountered: rss_ringoccs
                        \r\t\t%s\n
                        \r\tCould not determine occultation type.
                        \r\tCheck your input GEO file.
                    """ % fname
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

        geo_rho = self.rho_km_vals[crange]
        self.D_km_vals = self.D_km_vals[crange]
        self.rho_dot_kms_vals = self.rho_dot_kms_vals[crange]
        self.phi_rad_vals = self.phi_rad_vals[crange]
        self.B_rad_vals = self.B_rad_vals[crange]
        self.rho_km_vals = np.arange(np.min(geo_rho),
                                     np.max(geo_rho), dx_km_desired)

        if verbose:
            print("Interpolating Data...")

        interp = interpolate.interp1d(geo_rho, self.D_km_vals)
        self.D_km_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.rho_dot_kms_vals)
        self.rho_dot_kms_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rad_vals)
        self.phi_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.B_rad_vals)
        self.B_rad_vals = interp(self.rho_km_vals)
        del geo_rho, interp

 'history'
 'p_norm_vals'
 'phase_rad_vals'
 'phi_rl_rad_vals'
 'power_vals'
 'raw_tau_threshold_vals'
 'rev_info'
 'rho_corr_pole_km_vals'
 'rho_corr_timing_km_vals'
 't_oet_spm_vals'
 't_ret_spm_vals'
 't_set_spm_vals'

        r0       = self.rho_km_vals
        b        = self.B_rad_vals
        d        = self.D_km_vals
        phi      = self.phi_rad_vals
        rng      = [rho-10.0*res,rho+10.0*res]
        nstar    = np.min((r0-rho>=0).nonzero())
        r        = r0[nstar]
        phistar  = phi[nstar]
        F = special_functions.fresnel_scale(lambda_km, self.D_km_vals,
                                            self.phi_rad_vals, self.B_rad_vals)

        T_hat_vals = (1.0-1.0j)*np.exp(1j*psi_vals)/(2.0*F)
        p_norm_vals = np.abs(T_hat_vals)*np.abs(T_hat_vals)
        phase_rad_vals = -np.arctan2(np.imag(T_hat_vals),np.real(T_hat_vals))

        lambda_vals = np.zeros(np.size(self.rho_km_vals))+lambda_km
        self.p_norm_vals = p_norm_vals
        self.phase_rad_vals = phase_rad_vals
        self.f_sky_hz_vals = diffraction_correction.SPEED_OF_LIGHT_KM/lambda_vals
        self.nstar = nstar
        self.history        = "Delta Impulse DIffraction Model"

        recdata = DiffractionCorrection(self, res, rng=rng, wtype=wtype,
                                        fwd=fwd, norm=norm, verbose=verbose,
                                        bfac=bfac, psitype=psitype)

        self = recdata
        self.rho_star = rho


class DeltaImpulseDiffraction(object):
    def __init__(self, geo, lambda_km, res, rho, dx_km_desired=0.25,
                 occ=False, wtype='kb25', fwd=False, norm=True, bfac=True,
                 verbose=True, psitype='full', usefres=False):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.SquareWellFromGEO"
        check_type(verbose, bool, "verbose", fname)
        check_type(norm, bool, "norm", fname)
        check_type(bfac, bool, "bfac", fname)
        check_type(occ, str, "occ", fname)
        check_type(fwd, bool, "fwd", fname)
        check_type(geo, str, "geo", fname)
        check_type(wtype, str, "wtype", fname)
        check_type(psitype, str, "psitype", fname)
        check_type_and_convert(lambda_km, float, "lambda_km", fname)
        check_type_and_convert(res, float, "res", fname)
        check_type_and_convert(dx_km_desired, float, "dx_km_desired", fname)

        check_positive(res, "res", fname)
        check_positive(dx_km_desired, "dx_km_desired", fname)
        check_positive(lambda_km, "lambda_km", fname)

        data = CSV_tools.get_geo(geo, verbose=verbose)

        self.__retrieve_variables(data,verbose)
        self.__compute_variables(dx_km_desired,occ,verbose)
        self.__interpolate_variables(verbose)

        r0       = self.rho_km_vals
        b        = self.B_rad_vals
        d        = self.D_km_vals
        phi      = self.phi_rad_vals
        rng      = [rho-10.0*res,rho+10.0*res]
        nstar    = np.min((r0-rho>=0).nonzero())
        r        = r0[nstar]
        phistar  = phi[nstar]
        F        = fresnel_scale(lambda_km,d,phi,b)

        if (usefres == True):
            psi_vals = (np.pi/2.0)*((r0-r)/F)*((r0-r)/F)
        else:
            kD          = 2.0 * np.pi * d[nstar] / lambda_km
            dphi_s  = psi_factor(r0,r,b,phi)
            phi_s       = phi - dphi_s
            loop        = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s)) > 1.e-8):
                psi_d1  = psi_d1_phi(r,r0,d,b,phi_s,phi0,error_check=False)
                psi_d2  = psi_d2_phi(r,r0,d,b,phi_s,phi0,error_check=False)
                dphi_s  = -psi_d1 / psi_d2
                phi_s  += dphi_s
                loop   += 1
                if loop > 20:
                    break
                if verbose: print("Psi Iter: %d" % loop,end="\r")
            if verbose: print("Psi Iter: %d  dphi: %" % (loop,np.max(dphi_s)))
            psi_vals = kD * psi(r,r0,d,b,phi_s,phi0)

        T_hat_vals     = (1.0-1.0j)*np.exp(1j*psi_vals)/(2.0*F)
        p_norm_vals    = np.abs(T_hat_vals)*np.abs(T_hat_vals)
        phase_rad_vals = -np.arctan2(np.imag(T_hat_vals),np.real(T_hat_vals))

        lambda_vals         = np.zeros(np.size(self.rho_km_vals))+lambda_km
        self.p_norm_vals    = p_norm_vals
        self.phase_rad_vals = phase_rad_vals
        self.f_sky_hz_vals  = SPEED_OF_LIGHT_KM/lambda_vals
        self.nstar          = nstar
        self.history        = "Delta Impulse DIffraction Model"

        recdata = diffraction_correction(self, res, rng=rng, wtype=wtype,
                                         fwd=fwd, norm=norm, verbose=verbose,
                                         bfac=bfac, psitype=psitype)

        self = recdata
        self.rho_star = rho

    def __retrieve_variables(self,geo_dat,verbose):
        if verbose: print("Retrieving Variables...")
        geo_rho          = np.array(geo_dat.rho_km_vals)
        geo_D            = np.array(geo_dat.D_km_vals)
        geo_drho         = np.array(geo_dat.rho_dot_kms_vals)
        geo_phi          = np.array(geo_dat.phi_ora_deg_vals)
        geo_B            = np.array(geo_dat.B_deg_vals)
        rhotype,drhotype = check_real(geo_rho),check_real(geo_drho)
        phitype,Btype    = check_real(geo_phi),check_real(geo_B)
        Dtype            = check_real(geo_D)
        if not rhotype:
            raise TypeError("Bad GEO: rho_km_vals not real valued.")
        elif (np.min(geo_rho) < 0.0):
            raise ValueError("Bad GEO: rho_km_vals has negative values.")
        else: self.geo_rho  = geo_rho
        if not Dtype:
            raise TypeError("Bad GEO: D_km_vals not real valued.")
        elif (np.min(geo_D) < 0.0):
            raise ValueError("Bad GEO: D_km_vals has negative values.")
        else: self.geo_D    = geo_D
        if not drhotype:
            raise TypeError("Bad GEO: rho_dot_kms_vals not real valued.")
        else: self.geo_drho = geo_drho
        if not phitype:
            raise TypeError("Bad GEO: phi_deg_ora_vals not real valued.")
        elif (np.max(np.abs(geo_phi)) > 360.0):
            raise ValueError("Bad GEO: phi_deg_ora_vals > 360")
        else: self.geo_phi  = geo_phi
        if not Btype:
            raise TypeError("Bad GEO: B_deg_vals not real valued.")
        elif (np.max(np.abs(geo_B)) > 360.0):
            raise ValueError("Bad GEO: B_de2g_vals > 360")
        else: self.geo_B    = geo_B
        del rhotype,Dtype,drhotype,phitype,Btype

    def __compute_variables(self,dx_km_desired,occ,verbose):
        if verbose: print("Computing Variables...")
        geo_drho        = self.geo_drho

        if (occ == 'ingress'):  crange = (geo_drho < 0.0).nonzero()
        elif (occ == 'egress'): crange = (geo_drho > 0.0).nonzero()
        else:
            crange_e = (geo_drho > 0.0).nonzero()
            crange_i = (geo_drho < 0.0).nonzero()
            n_e      = np.size(crange_e)
            n_i      = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    "rho_dot_kms_vals has positive and negative values.\
                    This is likely a chord occultation. Set occ='ingress'\
                    to examine the ingress portion, and occ='egress'\
                    for the egress portion.")
            elif (n_e == 0) and (n_i == 0):
                raise ValueError("rho_dot_kms_vals is either empty or zero.")
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else: raise TypeError("Bad Input: GEO DATA")
        del n_e,n_i,crange_e,crange_i
        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else: raise ValueError("Bad occ input: Set 'egress' or 'ingress'")
            raise ValueError("Bad occ Input: '%s': %s" % (occ,mes))
        self.geo_rho        = self.geo_rho[crange]
        self.geo_D          = self.geo_D[crange]
        self.geo_drho       = self.geo_drho[crange]
        self.geo_phi        = self.geo_phi[crange]
        self.geo_B          = self.geo_B[crange]
        rmin                = np.min(geo_rho)
        rmax                = np.max(geo_rho)
        self.rho_km_vals    = np.arange(rmin,rmax,dx_km_desired)
        del rmin,rmax,geo_rho,geo_D,geo_drho,geo_phi,geo_B,rho_km_vals

    def __interpolate_variables(self,verbose):
        if verbose: print("Interpolating Data...")
        rho_km_vals           = self.rho_km_vals
        geo_rho               = self.geo_rho
        geo_drho              = self.geo_drho
        geo_D                 = self.geo_D
        geo_phi               = self.geo_phi
        geo_B                 = self.geo_B
        D_km_interp           = interpolate.interp1d(geo_rho,geo_D)
        self.D_km_vals        = D_km_interp(rho_km_vals)
        rho_dot_interp        = interpolate.interp1d(geo_rho,geo_drho)
        self.rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        phi_deg_interp        = interpolate.interp1d(geo_rho,geo_phi)
        phi_deg_vals          = phi_deg_interp(rho_km_vals)
        B_deg_interp          = interpolate.interp1d(geo_rho,geo_B)
        B_deg_vals            = B_deg_interp(rho_km_vals)
        self.phi_rad_vals     = np.deg2rad(phi_deg_vals)
        self.B_rad_vals       = np.deg2rad(B_deg_vals)
        del geo_rho,geo_drho,geo_D,geo_phi,geo_B,D_km_interp,rho_dot_interp
        del phi_deg_vals,phi_deg_interp,B_deg_interp,B_deg_vals
      