import numpy as np
from scipy import interpolate
from . import diffraction_correction, special_functions, window_functions
from rss_ringoccs.tools import CSV_tools, error_check
from rss_ringoccs.tools.history import write_history_dict, date_to_rev, rev_to_occ_info


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
                 verbose=True, psitype='fresnel', use_fresnel=False):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.SquareWellFromGEO"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(use_fresnel, bool, "use_fresnel", fname)
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
            print("\tRetrieving Variables...")

        try:
            # Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = np.array(data.rho_km_vals)
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = np.array(data.phi_ora_deg_vals)
            self.phi_rad_vals = np.deg2rad(self.phi_rad_vals)
            errmess = "B_rad_vals"
            self.B_rad_vals = np.array(data.B_deg_vals)
            self.B_rad_vals = np.deg2rad(self.B_rad_vals)
            errmess = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(data.t_oet_spm_vals)
            errmess = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(data.t_ret_spm_vals)
            errmess = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(data.t_set_spm_vals)
            errmess = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = np.array(data.phi_rl_deg_vals)
            self.phi_rl_rad_vals = np.deg2rad(self.phi_rl_rad_vals)
            errmess = "D_km_vals"
            self.D_km_vals = np.array(data.D_km_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(data.rho_dot_kms_vals)
        except (ValueError, TypeError, NameError, AttributeError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not convert %s into a numpy array.
                    \r\tCheck input GEO file.
                """ % (fname, errmess)
            )


        # Run error check on variables from GEO file.
        error_check.check_is_real(self.D_km_vals, "D_km_vals", fname)
        error_check.check_is_real(self.B_rad_vals, "B_rad_vals", fname)
        error_check.check_is_real(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_is_real(self.phi_rad_vals, "phi_rad_vals", fname)
        error_check.check_is_real(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_is_real(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_is_real(self.t_set_spm_vals, "t_set_spm_vals", fname)
        error_check.check_is_real(self.phi_rl_rad_vals,
                                  "phi_rl_rad_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals,
                                  "rho_dot_km_vals", fname)

        error_check.check_positive(self.D_km_vals, "D_km_vals", fname)
        error_check.check_positive(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_positive(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_positive(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_positive(self.t_set_spm_vals, "t_set_spm_vals", fname)

        error_check.check_two_pi(self.B_rad_vals, "B_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rad_vals, "phi_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rl_rad_vals, "phi_rl_rad_vals",
                                 fname, deg=False)

        if verbose:
            print("\tComputing Variables...")

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
            print("\tInterpolating Data...")

        interp = interpolate.interp1d(geo_rho, self.D_km_vals)
        self.D_km_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.rho_dot_kms_vals)
        self.rho_dot_kms_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rad_vals)
        self.phi_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rl_rad_vals)
        self.phi_rl_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.B_rad_vals)
        self.B_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_oet_spm_vals)
        self.t_oet_spm_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_ret_spm_vals)
        self.t_ret_spm_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_set_spm_vals)
        self.t_set_spm_vals = interp(self.rho_km_vals)
        del geo_rho, interp

        # Add 'fake' variables necessary for DiffractionCorrection class.
        self.raw_tau_threshold_vals = np.zeros(np.size(self.rho_km_vals))
        self.rho_corr_pole_km_vals = np.zeros(np.size(self.rho_km_vals))
        self.rho_corr_timing_km_vals = np.zeros(np.size(self.rho_km_vals))
        self.f_sky_hz_vals = np.zeros(np.size(self.rho_km_vals))
        self.f_sky_hz_vals += diffraction_correction.SPEED_OF_LIGHT_KM/lambda_km

        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data":     geo,
            "Wavelength":   lambda_km,
            "Resolution":   res,
            "Radius":       rho,
            "Width":        width
        }

        input_kwds = {
            "Sample Spacing":   dx_km_desired,
            "Occultation Type": occ,
            "Window Type":      wtype,
            "Forward Model":    fwd,
            "Normalization":    norm,
            "b-factor":         bfac,
            "verbose":          verbose,
            "Psi Type":         psitype,
            "Fresnel Model":    use_fresnel
        }

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        var = geo.split("/")[-1]

        try:
            var = var.split("_")
            band = '"%s"' % var[3][0]
            year = var[1]
            doy = var[2]
            dsn = "DSS-%s" % (var[3][1:])
            rev_num = date_to_rev(int(year), int(doy))
            occ_dir = rev_to_occ_info(rev_num)
            prof_dir = '"%s%"' % var[4]
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

        if use_fresnel:
            center = np.min((self.rho_km_vals >= rho).nonzero())
            cb2 = np.cos(self.B_rad_vals[center])
            cb2 *= cb2
            sb2 = np.sin(self.B_rad_vals[center])
            sb2 *= sb2
            sp2 = np.sin(self.phi_rad_vals[center])
            sp2 *= sp2
            F = np.sqrt(0.5*lambda_km*self.D_km_vals[center]*(1-cb2*sp2)/(sb2))
            T_hat = special_functions.square_well_diffraction(self.rho_km_vals,
                                                              rho-width/2.0,
                                                              rho+width/2.0, F)
            self.p_norm_vals = np.abs(T_hat)
            self.phase_rad_vals = np.arctan2(np.imag(T_hat), np.real(T_hat))
        else:
            self.p_norm_vals = np.zeros(np.size(self.rho_km_vals))+1.0
            rstart = np.min((self.rho_km_vals>=rho-width/2.0).nonzero())
            rfinsh = np.max((self.rho_km_vals<=rho+width/2.0).nonzero())
            self.p_norm_vals[rstart:rfinsh+1] = 0.0
            self.phase_rad_vals = np.zeros(np.size(self.rho_km_vals))
            rec = diffraction_correction.DiffractionCorrection(
                self, res, psitype=psitype, verbose=verbose, wtype=wtype,
                bfac=bfac
            )

            self.p_norm_vals = rec.power_vals
            self.phase_rad_vals = -rec.phase_vals

            crange = np.arange(rec.start, rec.start+rec.n_used, 1)
            self.B_rad_vals = self.B_rad_vals[crange]
            self.D_km_vals = self.D_km_vals[crange]
            self.f_sky_hz_vals = self.f_sky_hz_vals[crange]
            self.phi_rad_vals = self.phi_rad_vals[crange]
            self.phi_rl_rad_vals = self.phi_rl_rad_vals[crange]
            self.raw_tau_threshold_vals = self.raw_tau_threshold_vals[crange]
            self.rho_corr_pole_km_vals = self.rho_corr_pole_km_vals[crange]
            self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals[crange]
            self.rho_dot_kms_vals = self.rho_dot_kms_vals[crange]
            self.rho_km_vals = self.rho_km_vals[crange]
            self.t_oet_spm_vals = self.t_oet_spm_vals[crange]
            self.t_ret_spm_vals = self.t_ret_spm_vals[crange]
            self.t_set_spm_vals = self.t_set_spm_vals[crange]

        if verbose:
            print("\tData Extraction Complete.")


class DeltaImpulseDiffraction(object):
    def __init__(self, geo, lambda_km, res, rho, width, dx_km_desired=0.25,
                 occ="other", wtype='kb25', fwd=False, norm=True, bfac=True,
                 verbose=True, psitype='fresnel'):

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
                                                           "dx_km_desired",
                                                           fname)

        error_check.check_positive(res, "res", fname)
        error_check.check_positive(width, "width", fname)
        error_check.check_positive(dx_km_desired, "dx_km_desired", fname)
        error_check.check_positive(lambda_km, "lambda_km", fname)

        data = CSV_tools.get_geo(geo, verbose=verbose)
        occ = occ.replace(" ", "").replace("'", "").replace('"', "")
        occ = occ.lower()

        if verbose:
            print("\tRetrieving Variables...")

        try:
            # Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = np.array(data.rho_km_vals)
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = np.array(data.phi_ora_deg_vals)
            self.phi_rad_vals = np.deg2rad(self.phi_rad_vals)
            errmess = "B_rad_vals"
            self.B_rad_vals = np.array(data.B_deg_vals)
            self.B_rad_vals = np.deg2rad(self.B_rad_vals)
            errmess = "t_oet_spm_vals"
            self.t_oet_spm_vals = np.array(data.t_oet_spm_vals)
            errmess = "t_ret_spm_vals"
            self.t_ret_spm_vals = np.array(data.t_ret_spm_vals)
            errmess = "t_set_spm_vals"
            self.t_set_spm_vals = np.array(data.t_set_spm_vals)
            errmess = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = np.array(data.phi_rl_deg_vals)
            self.phi_rl_rad_vals = np.deg2rad(self.phi_rl_rad_vals)
            errmess = "D_km_vals"
            self.D_km_vals = np.array(data.D_km_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = np.array(data.rho_dot_kms_vals)
        except (ValueError, TypeError, NameError, AttributeError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tCould not convert %s into a numpy array.
                    \r\tCheck input GEO file.
                """ % (fname, errmess)
            )

        # Run error check on variables from GEO file.
        error_check.check_is_real(self.D_km_vals, "D_km_vals", fname)
        error_check.check_is_real(self.B_rad_vals, "B_rad_vals", fname)
        error_check.check_is_real(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_is_real(self.phi_rad_vals, "phi_rad_vals", fname)
        error_check.check_is_real(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_is_real(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_is_real(self.t_set_spm_vals, "t_set_spm_vals", fname)
        error_check.check_is_real(self.phi_rl_rad_vals,
                                  "phi_rl_rad_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals,
                                  "rho_dot_km_vals", fname)

        error_check.check_positive(self.D_km_vals, "D_km_vals", fname)
        error_check.check_positive(self.rho_km_vals, "rho_km_vals", fname)
        error_check.check_positive(self.t_oet_spm_vals, "t_oet_spm_vals", fname)
        error_check.check_positive(self.t_ret_spm_vals, "t_ret_spm_vals", fname)
        error_check.check_positive(self.t_set_spm_vals, "t_set_spm_vals", fname)

        error_check.check_two_pi(self.B_rad_vals, "B_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rad_vals, "phi_rad_vals",
                                 fname, deg=False)
        error_check.check_two_pi(self.phi_rl_rad_vals, "phi_rl_rad_vals",
                                 fname, deg=False)

        if verbose:
            print("\tComputing Variables...")

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
            print("\tInterpolating Data...")

        interp = interpolate.interp1d(geo_rho, self.D_km_vals)
        self.D_km_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.rho_dot_kms_vals)
        self.rho_dot_kms_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rad_vals)
        self.phi_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.phi_rl_rad_vals)
        self.phi_rl_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.B_rad_vals)
        self.B_rad_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_oet_spm_vals)
        self.t_oet_spm_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_ret_spm_vals)
        self.t_ret_spm_vals = interp(self.rho_km_vals)
        interp = interpolate.interp1d(geo_rho, self.t_set_spm_vals)
        self.t_set_spm_vals = interp(self.rho_km_vals)
        del geo_rho, interp

        # Add 'fake' variables necessary for DiffractionCorrection class.
        self.raw_tau_threshold_vals = np.zeros(np.size(self.rho_km_vals))
        self.rho_corr_pole_km_vals = np.zeros(np.size(self.rho_km_vals))
        self.rho_corr_timing_km_vals = np.zeros(np.size(self.rho_km_vals))
        self.f_sky_hz_vals = np.zeros(np.size(self.rho_km_vals))
        self.f_sky_hz_vals += diffraction_correction.SPEED_OF_LIGHT_KM/lambda_km

        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data":     geo,
            "Wavelength":   lambda_km,
            "Resolution":   res,
            "Radius":       rho,
            "Width":        width
        }

        input_kwds = {
            "Sample Spacing":   dx_km_desired,
            "Occultation Type": occ,
            "Window Type":      wtype,
            "Forward Model":    fwd,
            "Normalization":    norm,
            "b-factor":         bfac,
            "verbose":          verbose,
            "Psi Type":         psitype
        }

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        var = geo.split("/")[-1]

        try:
            var = var.split("_")
            band = '"%s"' % var[3][0]
            year = var[1]
            doy = var[2]
            dsn = "DSS-%s" % (var[3][1:])
            rev_num = date_to_rev(int(year), int(doy))
            occ_dir = rev_to_occ_info(rev_num)
            prof_dir = '"%s%"' % var[4]
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

        center = np.min((self.rho_km_vals >= rho).nonzero())
        cb2 = np.cos(self.B_rad_vals[center])
        cb2 *= cb2
        sb2 = np.sin(self.B_rad_vals[center])
        sb2 *= sb2
        sp2 = np.sin(self.phi_rad_vals[center])
        sp2 *= sp2
        F = np.sqrt(0.5*lambda_km*self.D_km_vals[center]*(1-cb2*sp2)/(sb2))
        T_hat = special_functions.square_well_diffraction(self.rho_km_vals,
                                                          rho-width/2.0,
                                                          rho+width/2.0, F)
        self.p_norm_vals = np.abs(T_hat)
        self.phase_rad_vals = np.arctan2(np.imag(T_hat), np.real(T_hat))

        if verbose:
            print("\tData Extraction Complete.")
   