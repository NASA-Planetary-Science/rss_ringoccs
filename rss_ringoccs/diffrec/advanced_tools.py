import numpy
from . import diffraction_correction, special_functions, window_functions
from rss_ringoccs.tools import CSV_tools, error_check, history


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
        wtype = error_check.check_wtype(wtype, fname)

        # Check that range and psitype are legal inputs.
        rng = error_check.check_range_input(rng, fname)
        psitype = error_check.check_psitype(psitype, fname)

        data = CSV_tools.ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
        rec = diffraction_correction.DiffractionCorrection(
            data, res, rng=rng, wtype=wtype, fwd=fwd, norm=norm, bfac=bfac,
            sigma=sigma, psitype=psitype, res_factor=res_factor, verbose=verbose
        )

        self.rho_km_vals = rec.rho_km_vals
        rmin = numpy.min(self.rho_km_vals)
        rmax = numpy.max(self.rho_km_vals)
        rstart = int(numpy.min((data.tau_rho-rmin>=0).nonzero()))
        rfin = int(numpy.max((rmax-data.tau_rho>=0).nonzero()))

        self.tau_rho = data.tau_rho[rstart:rfin+1]
        self.tau_power = data.power_vals[rstart:rfin+1]
        self.tau_phase = data.phase_vals[rstart:rfin+1]
        self.tau_tau = data.tau_vals[rstart:rfin+1]

        rmin = numpy.min(data.tau_rho)
        rmax = numpy.max(data.tau_rho)
        rstart = int(numpy.min((self.rho_km_vals-rmin>=0).nonzero()))
        rfin = int(numpy.max((rmax-self.rho_km_vals>=0).nonzero()))

        self.rho_km_vals = rec.rho_km_vals[rstart:rfin+1]
        self.power_vals = rec.power_vals[rstart:rfin+1]
        self.phase_vals = rec.phase_vals[rstart:rfin+1]
        self.tau_vals = rec.tau_vals[rstart:rfin+1]

        self.wtype = wtype
        self.res = res


class FindOptimalResolution(object):
    def __init__(self, geo, cal, dlp, tau, sres, dres, nres,
                 norm=True, bfac=True, sigma=2.e-13, psitype="fresnel4",
                 rng='all', wlst=['kbmd20'], res_factor=0.75, verbose=True):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.FindOptimalResolution"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(cal, str, "cal", fname)
        error_check.check_type(dlp, str, "dlp", fname)
        error_check.check_type(tau, str, "tau", fname)
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
        rng = error_check.check_range_input(rng, fname)

        if (not all(isinstance(x, str) for x in wlst)):
            raise TypeError(
                """
                \n\r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
                \r\twlst must be a list of strings.
                """ % fname
            )
        else:
            for i in range(len(wlst)):
                wlst[i] = error_check.check_wtype(wlst[i], fname)

        nwins = numpy.size(wlst)
        self.linf = numpy.zeros((nres, nwins))
        self.l2 = numpy.zeros((nres, nwins))
        self.ideal_res = numpy.zeros((nwins))
        eres = sres + (nres-1)*dres
        res = sres
        data = CSV_tools.ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
        rec = diffraction_correction.DiffractionCorrection(data, res, rng=rng,
                                                           wtype="kb35")
        start = int(numpy.min((data.tau_rho-numpy.min(rec.rho_km_vals)>=0).nonzero()))
        fin = int(numpy.max((numpy.max(rec.rho_km_vals)-data.tau_rho>=0).nonzero()))
        tau_power = data.power_vals[start:fin+1]
        for i in numpy.arange(nres):
            for j in range(nwins):
                wtype = wlst[j]
                recint = diffraction_correction.DiffractionCorrection(data, res, rng=rng, wtype=wtype)
                p_int = numpy.abs(recint.power_vals - tau_power)
                self.linf[i,j] = numpy.max(p_int)
                self.l2[i,j] = numpy.sqrt(numpy.sum(p_int*p_int)*recint.dx_km)

                if verbose:
                    print("%s %f %s %f %s %s"
                          % ('Res:',res,'Max:',eres,"WTYPE:",wtype))

            res += dres
        self.ideal_res = sres+dres*numpy.argmin(self.linf, axis=0)


class ModelFromGEO(object):
    def __init__(self, geo, lambda_km, res, rho, width=100, dx_km_desired=0.25,
                 wtype='kb25', norm=True, bfac=True, sigma=2.e-13, opacity=None,
                 verbose=True, psitype='fresnel', use_fresnel=False,
                 eccentricity=0.0, periapse=0.0, use_deprecate=False, rng="all",
                 model="squarewell", echo=False, rho_shift=0.0,
                 data_rho=None, data_pow=None, data_phase=None):

        # Check all input variables for errors.
        fname = "diffrec.advanced_tools.ModelFromGEO"
        error_check.check_type(verbose, bool, "verbose", fname)
        error_check.check_type(use_fresnel, bool, "use_fresnel", fname)
        error_check.check_type(norm, bool, "norm", fname)
        error_check.check_type(bfac, bool, "bfac", fname)
        error_check.check_type(geo, str, "geo", fname)
        error_check.check_type(echo, bool, "echo", fname)
        error_check.check_type(wtype, str, "wtype", fname)
        error_check.check_type(psitype, str, "psitype", fname)
        res = error_check.check_type_and_convert(res, float, "res", fname)
        width = error_check.check_type_and_convert(width, float, "width", fname)
        lambda_km = error_check.check_type_and_convert(lambda_km, float,
                                                       "lambda_km", fname)
        rho_shift = error_check.check_type_and_convert(rho_shift, float,
                                                       "rho_shift", fname)
        dx_km_desired = error_check.check_type_and_convert(dx_km_desired, float,
                                                           "dx_km_desired", fname)

        eccentricity = error_check.check_type_and_convert(eccentricity, float,
                                                          "eccentricity", fname)
        periapse = error_check.check_type_and_convert(periapse, float,
                                                      "periapse", fname)

        model = error_check.check_model(model, fname)

        error_check.check_positive(res, "res", fname)
        error_check.check_positive(width, "width", fname)
        error_check.check_positive(dx_km_desired, "dx_km_desired", fname)
        error_check.check_positive(lambda_km, "lambda_km", fname)

        if opacity:
            error_check.check_type(opacity, float, "opacity", fname)
            error_check.check_non_negative(opacity, "opacity", fname)
        else:
            opacity = 0.0

        # Get the data from the selected CSV files.
        data = CSV_tools.get_geo(geo, verbose=verbose,
                                 use_deprecate=use_deprecate)

        # If modeling from real data, check the inputs.
        if ((type(data_pow) == type(None)) and
            (not (type(data_rho) == type(None)))):
            raise TypeError(
                """
                \n\r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
                \r\tYou gave input for data_rho but not data_pow.
                \r\tPlease provide an input array for power.
                """ % (fname)
            )
        elif ((type(data_rho) == type(None)) and
              (not (type(data_pow) == type(None)))):
            raise TypeError(
                """
                \n\r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
                \r\tYou gave input for data_pow but not data_rho.
                \r\tPlease provide an input array for rho.
                """ % (fname)
            )
        elif not ((type(data_rho) == type(None)) and
                  (type(data_pow) == type(None))):
            error_check.check_is_real(data_rho, "data_rho", fname)
            error_check.check_is_real(data_pow, "data_pow", fname)
            error_check.check_non_negative(data_rho, "data_rho", fname)
            error_check.check_non_negative(data_pow, "data_pow", fname)
            error_check.check_lengths(data_rho, data_pow, "data_rho",
                                      "data_pow", fname)
            if not (type(data_phase) == type(None)):
                error_check.check_is_real(data_phase, "data_phase", fname)
                error_check.check_lengths(data_rho, data_phase, "data_rho",
                                          "data_phase", fname)
            else:
                data_phase = numpy.zeros(numpy.size(data_rho))
        else:
            pass

        if verbose:
            print("\tRetrieving Variables...")

        # Attempt conversion of data to numpy arrays.
        try:
            # Create dummy variable in case an error occurs.
            errmess = "rho_km_vals"
            self.rho_km_vals = numpy.array(data.rho_km_vals)
            errmess = "phi_ora_rad_vals"
            self.phi_rad_vals = numpy.array(data.phi_ora_deg_vals)
            self.phi_rad_vals = numpy.deg2rad(self.phi_rad_vals)
            errmess = "B_rad_vals"
            self.B_rad_vals = numpy.array(data.B_deg_vals)
            self.B_rad_vals = numpy.deg2rad(self.B_rad_vals)
            errmess = "t_oet_spm_vals"
            self.t_oet_spm_vals = numpy.array(data.t_oet_spm_vals)
            errmess = "t_ret_spm_vals"
            self.t_ret_spm_vals = numpy.array(data.t_ret_spm_vals)
            errmess = "t_set_spm_vals"
            self.t_set_spm_vals = numpy.array(data.t_set_spm_vals)
            errmess = "phi_rl_rad_vals"
            self.phi_rl_rad_vals = numpy.array(data.phi_rl_deg_vals)
            self.phi_rl_rad_vals = numpy.deg2rad(self.phi_rl_rad_vals)
            errmess = "D_km_vals"
            self.D_km_vals = numpy.array(data.D_km_vals)
            errmess = "rho_dot_kms_vals"
            self.rho_dot_kms_vals = numpy.array(data.rho_dot_kms_vals)
        except (ValueError, TypeError, NameError, AttributeError):
            raise TypeError(
                """
                \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
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
        error_check.check_is_real(self.phi_rl_rad_vals, "phi_rl_rad_vals", fname)
        error_check.check_is_real(self.rho_dot_kms_vals, "rho_dot_km_vals", fname)

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

        if res < 1.999999*dx_km_desired:
            raise ValueError(
                """
                \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
                \r\tResolution is less than twice the sample spacing.
                \r\t\tRequested Resolution (km):    %f
                \r\t\tSample Spacing (km):          %f\n
                \r\tChoose a resolution GREATER than %f km\n
                """ % (fname, res, dx_km_desired, 2.0*dx_km_desired)
            )

        if verbose:
            print("\tComputing Variables...")

        # Check that rho_km_vals is increasing and the rev isn't a chord occ.
        drho = [numpy.min(self.rho_dot_kms_vals), numpy.max(self.rho_dot_kms_vals)]
        dx_km = self.rho_km_vals[1] - self.rho_km_vals[0]

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
                \r\t\tSplit the input into two parts: Ingress and Engress
                \r\t\tand perform diffraction correction twice.
                """ % (fname)
            )
        elif ((drho[0] == 0.0) or (drho[1] == 0.0)):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tdrho/dt has zero valued elements.
                    \r\tYour input file is probably a chord occultation.
                    \r\tDiffraction Correction can only be performed for
                    \r\tone event at a time. That is, ingress or egress.\n
                    \r\tTO CORRECT THIS:
                    \r\t\tSplit the input into two parts: Ingress and Engress
                    \r\t\tand perform diffraction correction twice.
                    \r\t\tIgnore the region where drho/dt is close to zero.
                """ % (fname)
            )
        elif (dx_km > 0) and (drho[1] < 0):
            self.rho_dot_kms_vals = numpy.abs(self.rho_dot_kms_vals)
        elif (dx_km < 0) and (drho[0] > 0):
            raise ValueError(
                """
                    \r\tError Encountered:
                    \r\t\t%s\n
                    \r\trho_km_vals is decreasing yet rho_dot_kms_vals
                    \r\tis positiive. Check DLP class for errors.
                """ % (fname)
            )
        elif (dx_km < 0):
            self.rho_km_vals = self.rho_km_vals[::-1]
            self.phi_rad_vals = self.phi_rad_vals[::-1]
            self.B_rad_vals = self.B_rad_vals[::-1]
            self.D_km_vals = self.D_km_vals[::-1]
            self.rho_dot_kms_vals = numpy.abs(self.rho_dot_kms_vals[::-1])
        else:
            del drho, dx_km

        geo_rho = self.rho_km_vals
        self.rho_km_vals = numpy.arange(numpy.min(geo_rho),
                                        numpy.max(geo_rho), dx_km_desired)

        if verbose:
            print("\tInterpolating Data...")

        # Interpolate data to be spaced evenly in radius.
        r = self.rho_km_vals
        self.rho_dot_kms_vals = numpy.interp(r, geo_rho, self.rho_dot_kms_vals)
        self.phi_rl_rad_vals = numpy.interp(r, geo_rho, self.phi_rl_rad_vals)
        self.t_ret_spm_vals = numpy.interp(r, geo_rho, self.t_ret_spm_vals)
        self.t_oet_spm_vals = numpy.interp(r, geo_rho, self.t_oet_spm_vals)
        self.t_set_spm_vals = numpy.interp(r, geo_rho, self.t_set_spm_vals)
        self.phi_rad_vals = numpy.interp(r, geo_rho, self.phi_rad_vals)
        self.B_rad_vals = numpy.interp(r, geo_rho, self.B_rad_vals)
        self.D_km_vals = numpy.interp(r, geo_rho, self.D_km_vals)

        # Add 'fake' variables necessary for DiffractionCorrection class.
        n_pts = numpy.size(r)
        self.raw_tau_threshold_vals = numpy.zeros(n_pts)
        self.rho_corr_pole_km_vals = numpy.zeros(n_pts)
        self.rho_corr_timing_km_vals = numpy.zeros(n_pts)
        self.f_sky_hz_vals = (numpy.zeros(n_pts) +
            special_functions.frequency_to_wavelength(lambda_km))

        # Interpolate modeled data, if necessary.
        if not (type(data_rho) == type(None)):
            rstart = numpy.min((data_rho >= numpy.min(r)).nonzero())
            rfinsh = numpy.max((data_rho <= numpy.max(r)).nonzero())
            data_rho = data_rho[rstart:rfinsh]
            data_pow = data_pow[rstart:rfinsh]
            data_phase = data_phase[rstart:rfinsh]

            self.data_pow = numpy.interp(r, data_rho, data_pow)
        else:
            pass

        del n_pts, geo_rho, r

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
            "Window Type":      wtype,
            "Normalization":    norm,
            "b-factor":         bfac,
            "verbose":          verbose,
            "Psi Type":         psitype,
            "Fresnel Model":    use_fresnel
        }

        self.history = history.write_history_dict(input_vars,
                                                  input_kwds, __file__)
        var = geo.split("/")[-1]

        try:
            var = var.split("_")
            band = '"%s"' % var[3][0]
            year = var[1]
            doy = var[2]
            dsn = "DSS-%s" % (var[3][1:])
            rev_num = history.date_to_rev(int(year), int(doy))
            occ_dir = history.rev_to_occ_info(rev_num)
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

        center = numpy.min((self.rho_km_vals >= rho).nonzero())
        F = special_functions.fresnel_scale(
            lambda_km, self.D_km_vals[center],
            self.phi_rad_vals[center], self.B_rad_vals[center]
        ) + numpy.zeros(numpy.size(self.rho_km_vals))

        kD_vals = special_functions.wavelength_to_wavenumber(lambda_km)
        kD_vals *= self.D_km_vals

        # Compute the Normalized Equaivalent Width (See MTR86 Equation 20)
        norm_eq = window_functions.func_dict[wtype]["normeq"]

        # Compute the window width. (See MTR86 Equations 19, 32, and 33).
        self.w_km_vals, Prange = window_functions.window_width(
            res, norm_eq, self.f_sky_hz_vals, F, self.rho_dot_kms_vals,
            sigma, bfac=bfac
        )

        # From the requested range, extract array of the form [a, b]
        if (isinstance(rng, str)):
            self.rng = numpy.array(diffraction_correction.region_dict[rng])
        else:
            self.rng = numpy.array([numpy.min(rng), numpy.max(rng)])

        # Compute the smallest and largest allowed radii for reconstruction.
        crho = self.rho_km_vals[Prange]
        w = self.w_km_vals[Prange]
        rho_min = self.rho_km_vals[Prange]-self.w_km_vals[Prange]/2.0
        rho_max = self.rho_km_vals[Prange]+self.w_km_vals[Prange]/2.0

        wrange = Prange[numpy.where((rho_min >= numpy.min(crho)) &
                                 (rho_max <= numpy.max(crho)))]

        # Check that there is enough data for reconstruction.
        if (numpy.size(wrange) == 0):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tThe window width is too large.
                    \r\t\tMinimum Available Radius:         %f
                    \r\t\tMaximum Available Radius:         %f
                    \r\t\tMinimum Required Window Width:    %f
                    \r\t\tMaximum Required Window Width:    %f
                """ % (fname, numpy.min(crho), numpy.max(crho),
                       numpy.min(w), numpy.max(w))
            )
        elif (numpy.max(crho) < numpy.min(self.rng)):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tMinimum requested range is out of bounds.
                    \r\tYour Requested Minimum (km):    %f
                    \r\tYour Requested Maximum (km):    %f
                    \r\tMaximum Available Data (km):    %f
                """ % (fname, numpy.min(self.rng),
                       numpy.max(self.rng), numpy.max(crho))
            )
        elif (numpy.min(crho) > numpy.max(self.rng)):
            raise ValueError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tMaximum requested range iss out of bounds.
                    \r\tYour Requested Minimum (km): %f
                    \r\tYour Requested Maximum (km): %f
                    \r\tMinimum Available Data (km): %f\n
                """ % (fname, numpy.min(self.rng),
                       numpy.max(self.rng), numpy.min(crho))
            )
        else:
            rho_min = numpy.min(crho[wrange])
            rho_max = numpy.max(crho[wrange])

            wrange = wrange[numpy.where((crho[wrange] >= numpy.min(self.rng)) &
                                    (crho[wrange] <= numpy.max(self.rng)))]

        if (numpy.size(wrange) <= 1):
            raise IndexError(
                """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\tRequested range is beyond available data.
                \r\t\tMinimum Possible Radius: %f
                \r\t\tMaximum Possible Radius: %f
                \r\t\tMinimum Requested Range: %f
                \r\t\tMaximum Requested Range: %f
                """ % (fname, rho_min, rho_max,
                       numpy.min(self.rng), numpy.max(self.rng))
            )
        else:
            start = wrange[0]
            finish = wrange[-1]
            n_used = 1 + (finish - start)

        if verbose:
            print("\tComputing Forward Model...")

        if (model == "squarewell"):
            self.p_norm_actual_vals = numpy.zeros(numpy.size(self.rho_km_vals))
            self.p_norm_actual_vals += 1.0
            rstart = numpy.min((self.rho_km_vals>=rho-width/2.0).nonzero())
            rfinsh = numpy.max((self.rho_km_vals<=rho+width/2.0).nonzero())
            self.p_norm_actual_vals[rstart:rfinsh+1] = opacity

            if use_fresnel:
                center = numpy.min((self.rho_km_vals >= rho).nonzero())
                T_hat = special_functions.square_well_diffraction(
                    self.rho_km_vals, rho-width/2.0, rho+width/2.0, F[center]
                )
                if opacity:
                    T_hat += special_functions.inverse_square_well_diffraction(
                        self.rho_km_vals, rho-width/2.0,
                        rho+width/2.0, F[center]
                    ) * numpy.sqrt(opacity)

        elif (model == "rightstraightedge"):
            self.p_norm_actual_vals = numpy.zeros(numpy.size(self.rho_km_vals))
            self.p_norm_actual_vals += opacity
            rstart = numpy.min((self.rho_km_vals>=rho).nonzero())
            self.p_norm_actual_vals[rstart:-1] = 1.0

            if use_fresnel:
                T_hat = special_functions.right_straightedge(self.rho_km_vals,
                                                             rho, F[center])
                if opacity:
                    T_hat += special_functions.left_straightedge(
                        self.rho_km_vals, rho, F[center]
                    ) * numpy.sqrt(opacity)

        elif (model == "leftstraightedge"):
            self.p_norm_actual_vals = numpy.zeros(numpy.size(self.rho_km_vals))
            self.p_norm_actual_vals += opacity
            rfinsh = numpy.max((self.rho_km_vals<=rho).nonzero())
            self.p_norm_actual_vals[0:rfinsh] = 1.0

            if use_fresnel:
                center = numpy.max((self.rho_km_vals <= rho).nonzero())
                T_hat = special_functions.left_straightedge(self.rho_km_vals,
                                                            rho, F[center])

                if opacity:
                    T_hat += special_functions.right_straightedge(
                        self.rho_km_vals, rho, F[center]
                    ) * numpy.sqrt(opacity)

        elif (model == "deltaimpulse"):
            center = numpy.min((self.rho_km_vals >= rho).nonzero())
            self.p_norm_actual_vals = numpy.zeros(numpy.size(self.rho_km_vals))
            self.p_norm_actual_vals[center] = 1.0/dx_km_desired
            if use_fresnel:
                psi = (numpy.pi/2.0)*numpy.square((rho-self.rho_km_vals)/F[center])
            else:
                # Compute Newton-Raphson perturbation
                kD = kD_vals[center]
                r0 = self.rho_km_vals
                r = r0[center]
                phi0 = self.phi_rad_vals
                phi = phi0[center]+numpy.zeros(numpy.size(r0))
                B = self.B_rad_vals[center]
                D = self.D_km_vals[center]

                psi_d1 = special_functions.fresnel_dpsi_dphi(
                    kD, r, r0, phi, phi0, B, D
                )
                loop = 0
                while (numpy.max(numpy.abs(psi_d1)) > 1.0e-4):
                    psi_d1 = special_functions.fresnel_dpsi_dphi(
                        kD, r, r0, phi, phi0, B, D
                    )
                    psi_d2 = special_functions.fresnel_d2psi_dphi2(kD, r, r0,
                                                                   phi, phi0,
                                                                   B, D)

                    # Newton-Raphson
                    phi += -(psi_d1 / psi_d2)

                    # Add one to loop variable for each iteration
                    loop += 1
                    if (loop > 3):
                        break
                psi = special_functions.fresnel_psi(kD, r, r0, phi, phi0, B, D)

            T_hat = numpy.exp(1j*psi)*(0.5-0.5j)/F[center]
            self.r = self.rho_km_vals
            self.psi = psi
        else:
            use_fresnel = False
            self.p_norm_vals = self.data_pow
            self.p_norm_actual_vals = self.data_pow

        if ((not use_fresnel) and (not (model == "deltaimpulse"))):
            T_in = numpy.sqrt(self.p_norm_actual_vals.astype(complex))
            T_hat = special_functions.fresnel_transform(
                T_in, self.rho_km_vals, F, self.w_km_vals,
                start, n_used, wtype, norm, True, psitype,
                self.phi_rad_vals, kD_vals, self.B_rad_vals, self.D_km_vals,
                periapse, eccentricity
            )
        else:
            T_hat = T_hat[start:start+n_used+1]

        if verbose:
            print("\tForward Model Complete.")

        self.p_norm_vals = numpy.abs(T_hat)*numpy.abs(T_hat)
        self.phase_rad_vals = -numpy.arctan2(numpy.imag(T_hat),
                                             numpy.real(T_hat))

        crange = numpy.arange(n_used+1)+start
        self.B_rad_vals = self.B_rad_vals[crange]
        self.D_km_vals = self.D_km_vals[crange]
        self.f_sky_hz_vals = self.f_sky_hz_vals[crange]
        self.p_norm_actual_vals = self.p_norm_actual_vals[crange]
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
        self.w_km_vals = self.w_km_vals[crange]

        if echo:
            if verbose:
                print("\tComputing Echo Model...")

            F = F[crange]
            kD_vals = kD_vals[crange]
            n_shift = int(rho_shift/dx_km_desired)
            self.phase_rad_vals = -numpy.roll(self.phase_rad_vals, n_shift)


            # Compute the smallest and largest allowed radii for reconstruction.
            erho = self.rho_km_vals
            rho_min = self.rho_km_vals-self.w_km_vals/2.0
            rho_max = self.rho_km_vals+self.w_km_vals/2.0

            wrange = Prange[numpy.where((rho_min >= numpy.min(erho)) &
                                    (rho_max <= numpy.max(erho)))]

            start = wrange[0]
            finish = wrange[-1]
            n_used = 1 + (finish - start)

            T_in = (self.p_norm_actual_vals.astype(complex) *
                    numpy.exp(1.0j*self.phase_rad_vals))

            if not (type(data_phase) == type(None)):
                self.data_phase = numpy.interp(self.rho_km_vals, data_rho,
                                               data_phase)
                self.phase_rad_vals -= self.data_phase

            T_hat = special_functions.fresnel_transform(
                T_in, self.rho_km_vals, F, self.w_km_vals,
                start, n_used, wtype, norm, False, psitype,
                self.phi_rad_vals, kD_vals, self.B_rad_vals, self.D_km_vals,
                periapse, eccentricity
            )

            self.p_norm_vals = numpy.square(numpy.abs(T_hat))
            self.phase_rad_vals = numpy.arctan2(numpy.imag(T_hat), numpy.real(T_hat))

            crange = numpy.arange(n_used+1)+start
            self.B_rad_vals = self.B_rad_vals[crange]
            self.D_km_vals = self.D_km_vals[crange]
            self.f_sky_hz_vals = self.f_sky_hz_vals[crange]
            self.p_norm_actual_vals = self.p_norm_actual_vals[crange]
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
            self.w_km_vals = self.w_km_vals[crange]

            if verbose:
                print("\tEcho Model Complete.")

        if verbose:
            print("\tData Extraction Complete.")
