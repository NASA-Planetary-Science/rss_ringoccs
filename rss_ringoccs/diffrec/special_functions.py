import numpy as np
from scipy.special import erf, lambertw
from . import window_functions

# Declare constants for multiples of pi.
TWO_PI = 6.283185307179586476925287
ONE_PI = 3.141592653589793238462643
SQRT_PI_2 = 1.253314137315500251207883
RADS_PER_DEGS = 0.0174532925199432957692369

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
        Purpose:
            To smooth data with a Savitzky-Golay filter.
            This removes high frequency noise while
            maintaining many of the original features of
            the input data.
        Arguments:
            :y:
                The input "Noisy" data.
            :window_size (*int*):
                The length of the window.
                Must be an odd number.
            :order (*int*):
                The order of the polynomial used for filtering.
                Must be less then window_size - 1.
        Keywords:
            :deriv (*int*):
                The order of the derivative what will be computed.
        Output:
            :y_smooth:
                The data smoothed by the Savitzky-Golay filter.
                This returns the nth derivative if the deriv
                keyword has been set.
        Notes:
            The Savitzky-Golay is a type of low-pass filter,
            particularly suited for smoothing noisy data.
            The main idea behind this approach is to make for
            each point a least-square fit with a polynomial of
            high order over a odd-sized window centered at the point.
        References:
            #. A. Savitzky, M. J. E. Golay, Smoothing and
               Differentiation of Data by Simplified Least Squares
               Procedures. Analytical Chemistry, 1964, 36 (8),
               pp 1627-1639.
            #. Numerical Recipes 3rd Edition: The Art of
               Scientific Computing W.H. Press, S.A. Teukolsky,
               W.T. Vetterling, B.P. Flannery Cambridge University
               Press ISBN-13: 9780521880688
        Dependencies:
            #. Numpy
    """
    try:
        y = np.array(y)
    except:
        raise TypeError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffrec Subpackage\n"
            "\tspecial_functions.savitzky_golay\n"
            "\tinput variable should be a numpy array."
            "\tSyntax:\n"
            "\t\tysmooth = savitzky_golay(y, Window_Size, Poly_Order)"
        )
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except (ValueError, TypeError):
        raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffrec Subpackage\n"
            "\tspecial_functions.savitzky_golay:\n"
            "\t\twindow_size must be an odd integer.\n"
        )

    if (window_size % 2 != 1) or (window_size < 1):
        raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffrec Subpackage\n"
            "\tspecial_functions.savitzky_golar:\n"
            "\t\twindow_size must be an odd integer.\n"
        )
    elif (window_size < order + 2):
        raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffrec Subpackage\n"
            "\tspecial_functions.savitzky_golar:\n"
            "\t\twindow_size is too small for the\n"
            "\t\trequested polynomial order.\n"
        )
    else:
        pass

    half_window = (window_size - 1) // 2

    # precompute coefficients
    b = np.zeros((window_size, order+1))
    b[..., 0] = 1
    for k in range(half_window):
        n0 = (half_window) - k
        m = n0
        n = -n0
        for j in range(1, order+1):
            b[k, j] = n
            b[window_size-1-k, j] = m
            n *= -n0
            m *= n0

    b = np.mat(b)

    m = 1
    for i in range(deriv):
        m *= rate*(deriv-i)

    m *= np.linalg.pinv(b).A[deriv]

    # Pad the endpoints with values from the signal.
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))

    return np.convolve(m[::-1], y, mode='valid')

def fresnel_transform(rho_km_vals, phi_rad_vals, F_km_vals, B_rad_vals,
                      d_km_vals, T_in, lambda_km_vals, w_km_vals, dx_km,
                      wtype="kb25", norm=True, fft=False, fwd=False,
                      verbose=True, psitype='full'):
    """
        Purpose:
            Computes the Fresnel Transform (Forward or Inverse)
            of a given data set.
        Arguments:
            rho_km_vals:
                Real Numpy Array
                Ring radius, in kilometers.
            phi_rad_vals:
                Real Numpy Array
                Ring azimuth angle, in radians.
            F_km_vals:
                Real Numpy Array
                Fresnel scale, in kilometers.
            B_rad_vals:
                Real Numpy Array
                Ring opening angle, in radians.
            lambda_sky_km_vals:
                Real Numpy Array
                Wavelength of recieved signal, in kilometers.
            D_km_vals:
                Real Numpy Array
                Spacecraft-RIP distance, in kilometers.
            dx_km:
                Float
                Sample spacing, in kilometers.
            T_in:
                Complex Numpy Array
                Complex data set.
            w_vals:
                Real Numpy Array
                Window width, in kilometers.
        Keywords:
            wtype:
                String
                Window used in reconstruction.
            norm:
                Boolean
                Boolean for determining if the transform
                will be normalized by the window width.
            fft:
                Boolean
                Whether or not to use FFT's.
            fwd:
                Boolean
                Used for determining if the Forward or
                Inverse Fresnel Transform will be computed.
                By default, the inverse transform is computed.
            verbose:
                Boolean
                Determines if progress reports print out.
            psitype:
                String
                Chosen approximation for the psi variable.
                Allowed strings are:
                    'full'      No Approximation is applied.
                    'MTR2'      Second Order Series from MTR86.
                    'MTR3'      Third Order Series from MTR86.
                    'MTR4'      Fourth Order Series from MTR86.
                    'Fresnel'   Standard Fresnel approximation.
        Output:
            T_out:
                The Fresnel transform of T_in. Forward transform
                if fwd=True was set, inverse transform otherwise.
        History:
            Translated from IDL: RJM - 2018/05/14 5:06 P.M.
            Updated: RJM - 2018/09/25 10:27 P.M.
    """
    # Retrieve requested window type, starting point, and sample spacing.
    start, n_used = window_functions.get_range_actual(rho_km_vals, [0.1, 1.e6],
                                                      w_km_vals)

    # Compute product of wavenumber and RIP distance.
    kD_vals = TWO_PI * d_km_vals / lambda_km_vals

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
    fw = window_functions.func_dict[wtype]["func"]

    # Define normalization function and verbose message.
    nrm = window_functions.normalize
    mes = "\t\tPt: %d  Tot: %d  Width: %d  Psi Iters: %d"

    # Set inverse function to FFT or Integration.
    finv = fresnel_inverse

    # If forward transform, adjust starting point by half a window.
    if fwd:
        w_max = np.max(w_km_vals[start:start + n_used])
        nw_fwd = int(np.ceil(w_max / (2.0 * dx_km)))
        start_fwd = int(start + nw_fwd)
        n_fwd = int(n_used - 2 * nw_fwd)
        start = start_fwd
        n_used = n_fwd
    else:
        pass
    
    # Create empty array for reconstruction / forward transform.
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
        r0 = rho_km_vals[center]
        r = rho_km_vals[crange]
        F2 = F_km_vals*F_km_vals
        x = (r-r0)
        psi_vals = (ONE_PI / 2.0) * x * x
        loop = 0
        for i in np.arange(n_used):
            # Current point being computed.
            center = start+i

            # Window width and Frensel scale for current point.
            w = w_km_vals[center]
            F = F_km_vals[center]

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
                r0 = rho_km_vals[center]

                # Compute psi for with stationary phase value
                x = r-r0
                psi_vals = (ONE_PI / 2.0) * x * x / F2[center]
            else:
                crange += 1
                psi_vals = (ONE_PI / 2.0) * x * x / F2[center]

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
        if verbose:
            print("\n", end="\r")
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
                print(mes % (i, n_used-1, nw, loop), end="\r")
        if verbose:
            print("\n", end="\r")
    return T_out

def compute_norm_eq(w_func, error_check=True):
    """
        Purpose:
            Compute normalized equivalenth width of a given function.
        Variables:
            :w_func:
                A numpy array. Function with which to compute
                the normalized equivalent width.
        Outputs:
            :normeq:
                The normalized equivalent width of w_func.
        Dependencies:
            #. numpy
        Notes:
            The normalized equivalent width is effectively computed
            using Riemann sums to approximate integrals. Therefore
            large dx values (Spacing between points in w_func)
            will result in an inaccurate normeq. One should keep
            this in mind during calculations.
        References:
            #. Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
               Profiling Saturn's rings by radio occultation,
               Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
               https://doi.org/10.1016/0019-1035(86)90078-3
        Examples:
            Compute the Kaiser-Bessel 2.5 window of width 30
            and spacing 0.1, and then use compute_norm_eq
            to compute the normalized equivalent width:
                In [1]: from rss_ringoccs import diffrec as dc
                In [2]: w = dc.window_functions.kb25(30, 0.1)
                In [3]: normeq = dc.special_functions.compute_norm_eq(w)
                In [4]: print(normeq)
                1.6573619266424229
            In contrast, the actual value is 1.6519208.
            Compute the normalized equivalent width for the squared
            cosine window of width 10 and spacing 0.25.
                In [1]: from rss_ringoccs import diffrec as dc
                In [2]: w = dc.window_functions.coss(10, 0.25)
                In [3]: normeq = dc.special_functions.compute_norm_eq(w)
                In [4]: print(normeq)
                1.5375000000000003
            The normalized equivalent width of the squared cosine
            function can be computed exactly using standard methods
            from a calculus course. It's value is exactly 1.5
            If we use a smaller dx when computing w, we get a better
            approximation. Use width 10 and spacing 0.001.
                In [1]: from rss_ringoccs import diffrec as dc
                In [2]: w = dc.window_functions.coss(10, 0.001)
                In [3]: normeq = dc.special_functions.compute_norm_eq(w)
                In [4]: print(normeq)
                1.50015
        History:
            Created: RJM - 2018/05/16 3:54 P.M.
    """
    if error_check:
        try:
            w_func = np.array(w_func)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.compute_norm_eq:\n"
                "\t\tFirst input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(w_func))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.compute_norm_eq:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point numbers.\n"
            )
        elif (np.min(w_func) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.compute_norm_eq:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass
    nw = np.size(w_func)
    tot = np.sum(w_func)
    normeq = nw*(np.sum(w_func*w_func)) / (tot*tot)

    return normeq

def fresnel_scale(Lambda, d, phi, b, deg=False, error_check=True):
    """
        Purpose:
            Compute the Fresnel Scale from lambda, D, Phi, and B.
        Variables:
            :Lambda (*array* or *float*):
                Wavelength of the incoming signal.
            :d (*array* or *float*):
                RIP-Spacecraft Distance.
            :phi (*array* or *float*):
                Ring azimuth angle.
            :b (*array* or *float*):
                Ring opening angle.
        Keywords:
            :deg (*bool*):
                Set True if phi/b are in degrees.
                Default is radians.
        Output:
            :FRES:
                The Fresnel scale.
        Note:
            Lambda and d must be in the same units.
            The output (Fresnel scale) will have the same units
            as lambda and d. In addition, b and phi must also
            have the same units. If b and phi are in degrees,
            make sure to set deg=True. Default is radians.
        History:
            Translated from IDL: RJM - 2018/04/15 12:36 P.M.
    """
    if error_check:
        try:
            Lambda = np.array(Lambda)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tFirst input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(Lambda))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point numbers.\n"
            )
        elif (np.min(Lambda) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass

        try:
            d = np.array(d)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tSecond input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(d))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tSecond input must be an array\n"
                "\t\tof floating point numbers.\n"
            )
        elif (np.min(d) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tSecond input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass

        try:
            phi = np.array(phi)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tThird input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(phi))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tThird input must be an array\n"
                "\t\tof floating point numbers.\n"
            )
        else:
            pass

        try:
            b = np.array(b)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tSecond input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(b))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.fresnel_scale:\n"
                "\t\tSecond input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass
    else:
        pass

    if deg:
        cb = np.cos(b * RADS_PER_DEGS)
        sb = np.sin(b * RADS_PER_DEGS)
        sp = np.sin(phi * RADS_PER_DEGS)
    else:
        cb = np.cos(b)
        sb = np.sin(b)
        sp = np.sin(phi)

    fres = np.sqrt(0.5 * Lambda * d * (1 - (cb*cb) * (sp*sp)) / (sb*sb))
    return fres

def psi_d1_phi(r, r0, d, b, phi, phi0, error_check=True):
    """
        Purpose:
            Calculate dpsi/dphi from geometric variables.
        Variables:
            :r:
                Numpy Array. Ring radius variable, in kilometers.
            :r0:
                Float. Ring intercept point, in kilometers.
            :d:
                Numpy Array. RIP-Spacecraft distance in kilometers.
            :b:
                Numpy Array. The ring opening angle, in radians.
            :phi:
                Numpy Array. The ring azimuth angle or r, in radians.
            :phi0:
                Float or Numpy Array. The ring azimuth angle
                for r0, in radians.
        History:
            Translated from IDL: RJM - 2018/05/15 7:06 P.M.
    """
    if error_check:
        try:
            r = np.array(r)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(r))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(r) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass
        
        if (not isinstance(r0, float)):
            try:
                r0 = float(r0)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffrec Subpackage\n"
                    "\tspecial_functions.psi_d1_phi:\n"
                    "\t\tSecond input must be a positive\n"
                    "\t\tfloating point number.\n"
                )
        else:
            pass
        
        if (np.min(r0) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tSecond input must be a positive\n"
                "\t\tfloating point number.\n"
            )
        else:
            pass

        try:
            d = np.array(d)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(d))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(d) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass

        try:
            b = np.array(b)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFourth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(b))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFourth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        try:
            phi = np.array(phi)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFifth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(phi))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFifth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        if (not isinstance(phi0, float)):
            try:
                phi0 = float(phi0)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffrec Subpackage\n"
                    "\tspecial_functions.psi_d1_phi:\n"
                    "\t\tSecond input must be a positive\n"
                    "\t\tfloating point number.\n"
                )
        else:
            pass
    else:
        pass

    cb = np.cos(b)
    sp = np.sin(phi)
    cp = np.cos(phi)
    sp0 = np.sin(phi0)
    cp0 = np.cos(phi0)
    xi = (cb / d) * (r0*cp0 - r*cp)
    eta = (r0*r0 + r*r - 2.0*r*r0*(sp*sp0 + cp*cp0)) / (d*d)
    v1 = r * cb * sp / d
    v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d*d)
    psi_d1_phi_vals = (2.0*v1 + v2) / (2.0 * np.sqrt(1.0 + 2.0*xi + eta)) - v1

    return psi_d1_phi_vals

def psi_d2_phi(r, r0, d, b, phi, phi0, error_check=True):
    """
        Purpose:
            Calculate dpsi/dphi from geometric variables.
        Variables:
            :r:
                Numpy Array. Ring radius variable, in kilometers.
            :r0:
                Float. Ring intercept point, in kilometers.
            :d:
                Numpy Array. RIP-Spacecraft distance in kilometers.
            :b:
                Numpy Array. The ring opening angle, in radians.
            :phi:
                Numpy Array. The ring azimuth angle or r, in radians.
            :phi0:
                Float or Numpy Array. The ring azimuth
                angle for r0, in radians.
    """
    if error_check:
        try:
            r = np.array(r)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFirst input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(r))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(r) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass
        
        if (not isinstance(r0, float)):
            try:
                r0 = float(r0)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffrec Subpackage\n"
                    "\tspecial_functions.psi_d2_phi:\n"
                    "\t\tSecond input must be a positive\n"
                    "\t\tfloating point number.\n"
                )
        else:
            pass
        
        if (np.min(r0) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tSecond input must be a positive\n"
                "\t\tfloating point number.\n"
            )
        else:
            pass

        try:
            r = np.array(d)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tThird input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(d))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(d) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass

        try:
            r = np.array(b)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFourth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(b))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFourth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        try:
            r = np.array(phi)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFifth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(phi))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d2_phi:\n"
                "\t\tFifth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        if (not isinstance(phi0, float)):
            try:
                phi0 = float(phi0)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffrec Subpackage\n"
                    "\tspecial_functions.psi_d2_phi:\n"
                    "\t\tSecond input must be a positive\n"
                    "\t\tfloating point number.\n"
                )
        else:
            pass
    else:
        pass

    cb = np.cos(b)
    sp = np.sin(phi)
    cp = np.cos(phi)
    sp0 = np.sin(phi0)
    cp0 = np.cos(phi0)
    xi = (cb / d) * (r0*cp0 - r*cp)
    eta = ((r0**2) + (r**2) - 2.0*r*r0*(sp*sp0 + cp*cp0)) / (d**2)
    v1 = r * cb * cp / d
    v2 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / (d**2)
    v3 = r * cb * sp / d
    v4 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d**2)
    dphia = (2.0*v1 + v2)/(2.0 * np.sqrt(1.0 + 2.0*xi + eta))
    dphib = v1 + ((2.0*v3 + v4)**2)/(4.0 * (np.sqrt(1.0 + 2.0*xi + eta)**3))
    psi_d2_phi_vals = dphia - dphib

    return psi_d2_phi_vals

def fresnel_inverse(T_hat, ker, dx, f_scale):
    """
        Purpose:
            Compute the approximation Fresnel
            Inverse (MTR86 Equation 15)
        Arguments:
            :T_hat:
                Complex Numpy Array.
                The complex transmittance of the
                normalized diffraction data.
            :ker:
                Complex Numpy Array. The Fresnel Kernel.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
            :f_scale:
                Real Numpy Array.
                The Fresnel Scale as a function of
                ring radius (km).
        Outputs:
            :T (*complex*):
                The fresnel inversion about the center
                of the Fresnel Kernel.
        Dependencies:
            #. numpy
    """
    T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
    return T

def psi_func(kD, r, r0, phi, phi0, B, D, error_check=True):
    """
        Purpose:
            Calculate psi from geometry variables.
        Arguments:
            :r:
                Numpy Array. Ring radius variable, in kilometers.
            :r0 (*float*):
                Ring intercept point, in kilometers.
            :d:
                Numpy Array. RIP-Spacecraft distance in kilometers.
            :b:
                Numpy Array. The ring opening angle, in radians.
            :phi:
                Numpy Array.
                The ring azimuth angle or r, in radians.
            :phi0:
                Float or Numpy Array.
                The ring azimuth angle for r0, in radians.
        Keywords:
            :verbose (*bool*):
                Boolean for printing out information to
                the command line.
    """
    if error_check:
        try:
            r = np.array(r)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(r))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(r) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass
        
        try:
            r0 = np.array(r0)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tSecond input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(r0))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(r0) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFirst input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass
        
        try:
            D = np.array(D)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(D))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof floating point number.\n"
            )
        elif (np.min(D) < 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tThird input must be an array\n"
                "\t\tof positive numbers.\n"
            )
        else:
            pass

        try:
            B = np.array(B)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFourth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(B))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFourth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        try:
            phi = np.array(phi)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFifth input could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        if (not np.all(np.isreal(phi))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffrec Subpackage\n"
                "\tspecial_functions.psi_d1_phi:\n"
                "\t\tFifth input must be an array\n"
                "\t\tof floating point number.\n"
            )
        else:
            pass

        if (not isinstance(phi0, float)):
            try:
                phi0 = float(phi0)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffrec Subpackage\n"
                    "\tspecial_functions.psi_d1_phi:\n"
                    "\t\tSecond input must be a positive\n"
                    "\t\tfloating point number.\n"
                )
        else:
            pass
    else:
        pass

    # Compute Xi variable (MTR86 Equation 4b).
    xi = (np.cos(B)/D) * (r * np.cos(phi) - r0 * np.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*np.cos(phi-phi0)) / (D*D)
    psi_vals = kD * (np.sqrt(1.0+eta-2.0*xi) - (1.0-xi))
    return psi_vals

def resolution_inverse(x, error_check=True):
    """
        Purpose:
            Compute the inverse of y = x/(exp(-x)+x-1)
        Variables:
            :x:
                A real or complex number, or numpy array.
        Outputs:
            :f:
                The inverse of x/(exp(-x)+x-1)
        Dependencies:
            #. numpy
            #. scipy.special
        Method:
            The inverse of x/(exp(-x)+x-1) is computed using the
            LambertW function. This function is the inverse of
            y = x * exp(x). This is computed using the scipy.special
            subpackage using their lambertw function.
        Warnings:
            #. The real part of the argument must be greater than 1.
            #. The scipy.special lambertw function is slightly
               inaccurate when it's argument is near -1/e. This
               argument is z = exp(x/(1-x)) * x/(1-x)
        References:
            [1] http://mathworld.wolfram.com/LambertW-Function.html
            [2] https://en.wikipedia.org/wiki/Lambert_W_function
        Examples:
            Plot the function on the interval (1,2)
                In [1]: import rss_ringoccs.diffcorr.special_functions as sf
                In [2]: import numpy as np
                In [3]: x = np.array(range(0,1001))*0.001+1.01
                In [4]: y = sf.resolution_inverse(x)
                In [5]: import matplotlib.pyplot as plt
                In [6]: plt.show(plt.plot(x,y))
                (Beautiful plots appear here)
    """
    if error_check:
        y = x
        try:
            x = np.array(x)
            if (not np.all(np.isreal(x))) and (not np.all(np.iscomplex(x))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\tspecial_function.resolution_inverse:\n"
                    "\t\tInput must be a real or complex valued numpy array.\n"
                    "\t\tThe elements of your array have type: %s"
                    % (x.dtype)
                )
            else:
                del y
        except (TypeError, ValueError) as errmes:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.resolution_inverse:\n"
                "\t\tInput must be a real or complex valued numpy array.\n"
                "\t\tYour input has type: %s\n"
                "\tOriginal Error Mesage: %s\n"
                % (type(y).__name__, errmes)
            )
        if (np.min(np.real(x)) <= 1.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.resolution_inverse:\n"
                "\t\tThis function is only defined for inputs\n"
                "\t\twhose real part is greater than 1.\n"
                "\t\tYour input has real minimum: %f\n"
                % (np.min(np.real(x)))
            )
        else:
            pass

    P1 = x/(1.0-x)
    P2 = P1*np.exp(P1)
    f = lambertw(P2)-P1

    if np.all(np.isreal(x)):
        f = np.real(f)

    return f

def fresnel_cos(x, error_check=True):
    """
        Purpose:
            Compute the Fresnel cosine function.
        Variables:
            :x:
                A real or complex number, or numpy array.
        Outputs:
            :f_cos:
                The fresnel cosine integral of x.
        Dependences:
            #. numpy
            #. scipy
        Notes:
            [1] The Fresnel Cosine integral is the solution to the
                equation dy/dx = cos(pi/2 * x^2), y(0) = 0. In other
                words, y = integral (t=0 to x) cos(pi/2 * t^2) dt
            [2] The Fresnel Cosine and Sine integrals are computed by
                using the scipy.special Error Function. The Error
                Function, usually denoted Erf(x), is the solution to
                dy/dx = (2/sqrt(pi)) * exp(-x^2), y(0) = 0. That is:
                y = 2/sqrt(pi) * integral (t=0 to x) exp(-t^2)dt.
                Using Euler's Formula for exponentials allows one
                to use this to solve for the Fresnel Cosine integral.
            [3] The Fresnel Cosine integral is used for the solution
                of diffraction through a square well. Because of this
                it is useful for forward modeling problems in 
                radiative transfer and diffraction.
        References:
            [1] https://en.wikipedia.org/wiki/Fresnel_integral
            [2] https://en.wikipedia.org/wiki/Error_function
            [3] http://mathworld.wolfram.com/FresnelIntegrals.html
            [4] http://mathworld.wolfram.com/Erf.html
        Examples:
            Compute and plot the Fresnel Cosine integral.
                In [1]: import rss_ringoccs.diffcorr.special_functions as sf
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = sf.fresnel_cos(x)
                In [6]: plt.show(plt.plot(x,y))
        History:
            Translated from IDL:     RJM - 2018/05/14 2:13 P.M.
            Allow complex arguments: RJM - 2018/05/16 1:25 P.M.
            Updated and added error check: RJM - 2018/09/19 7:00 P.M.
    """
    if error_check:
        y = x
        try:
            x = np.array(x)
            if (not np.all(np.isreal(x))) and (not np.all(np.iscomplex(x))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\tspecial_functions.fresnel_cos:\n"
                    "\t\tInput must be a real or complex valued numpy array.\n"
                    "\t\tThe elements of your array have type: %s"
                    % (x.dtype)
                )
            else:
                del y
        except (TypeError, ValueError) as errmes:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.fresnel_cos:\n"
                "\t\tInput must be a real or complex valued numpy array.\n"
                "\t\tYour input has type: %s\n"
                "\tOriginal Error Mesage: %s\n"
                % (type(y).__name__, errmes)
            )
    else:
        pass

    f_cos = ((0.25-0.25j)*erf((1.0+1.0j)*x*SQRT_PI_2)+
             (0.25+0.25j)*erf((1.0-1.0j)*x*SQRT_PI_2))

    if (np.isreal(x).all()):
        f_cos = np.real(f_cos)

    return f_cos

def fresnel_sin(x, error_check=True):
    """
        Purpose:
            Compute the Fresnel sine function.
        Variables:
            :x:
                A real or complex argument, or numpy array.
        Outputs:
            :f_sin:
                The fresnel sine integral of x.
        Dependences:
            #. numpy
            #. scipy
        Notes:
            [1] The Fresnel sine integral is the solution to the
                equation dy/dx = sin(pi/2 * x^2), y(0) = 0. In other
                words, y = integral (t=0 to x) sin(pi/2 * t^2) dt
            [2] The Fresnel Cossine and Sine integrals are computed
                by using the scipy.special Error Function. The Error
                Function, usually denoted Erf(x), is the solution to
                dy/dx = (2/sqrt(pi)) * exp(-x^2), y(0) = 0. That is:
                y = 2/sqrt(pi) * integral (t=0 to x) exp(-t^2)dt.
                Using Euler's Formula for exponentials allows one
                to use this to solve for the Fresnel Sine integral.
            [3] The Fresnel sine integral is used for the solution
                of diffraction through a square well. Because of this
                is is useful for forward modeling problems in 
                radiative transfer and diffraction.
        References:
            [1] https://en.wikipedia.org/wiki/Fresnel_integral
            [2] https://en.wikipedia.org/wiki/Error_function
            [3] http://mathworld.wolfram.com/FresnelIntegrals.html
            [4] http://mathworld.wolfram.com/Erf.html
        Examples:
            Compute and plot the Fresnel Sine integral.
                In [1]: import rss_ringoccs.diffcorr.special_functions as sf
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = sf.fresnel_sin(x)
                In [6]: plt.show(plt.plot(x,y))
    """
    if error_check:
        y = x
        try:
            x = np.array(x)
            if (not np.all(np.isreal(x))) and (not np.all(np.iscomplex(x))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\tspecial_function.fresnel_cos:\n"
                    "\t\tInput must be a real or complex valued numpy array.\n"
                    "\t\tThe elements of your array have type: %s"
                    % (x.dtype)
                )
            else:
                del y
        except (TypeError, ValueError) as errmes:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.fresnel_cos:\n"
                "\t\tInput must be a real or complex valued numpy array.\n"
                "\t\tYour input has type: %s\n"
                "\tOriginal Error Mesage: %s\n"
                % (type(y).__name__, errmes)
            )
    else:
        pass

    f_sin = ((0.25+0.25j)*erf((1.0+1.0j)*x*SQRT_PI_2)+
             (0.25-0.25j)*erf((1.0-1.0j)*x*SQRT_PI_2))

    if (np.isreal(x).all()):
        f_sin = np.real(f_sin)

    return f_sin

def single_slit_diffraction(x, z, a):
    """
        Purpose:
            Compute diffraction through a single slit for the
            variable x with a distance z from the slit and
            slit parameter a. This assume Fraunhofer diffraction.
        Variables:
            :x:
                A real or complex argument, or numpy array.
            :z (*float*):
                Float
                The perpendicular distance from the slit plane to
                the observer.
            :a (*float*):
                The slit parameter. This is a unitless paramter
                defined as the ratio between the slit width and
                the wavelength of the incoming signal.
        Outputs:
            :f:
                Single slit diffraction pattern.
        Dependences:
            [1] numpy
    """
    y = x
    try:
        x = np.array(x)
        if (not np.all(np.isreal(x))) and (not np.all(np.iscomplex(x))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.single_slit_diffraction_solve:\n"
                "\t\tInput must be a real or complex valued numpy array.\n"
                "\t\tThe elements of your array have type: %s"
                % (x.dtype)
            )
        else:
            del y
    except (TypeError, ValueError) as errmes:
        raise TypeError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.single_slit_diffraction_solve:\n"
            "\t\tFirst input must be a real or complex valued numpy array.\n"
            "\t\tYour input has type: %s\n"
            "\tOriginal Error Mesage: %s\n"
            % (type(y).__name__, errmes)
        )
    if not (isinstance(z, float)):
        try:
            z = float(z)
        except ValueError:
            raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.single_slit_diffraction_solve:\n"
            "\t\tSecond input must be a floating point number.\n"
            "\t\tYour input has type: %s\n"
            % (type(z).__name__)
            )
    if not (isinstance(a, float)):
        try:
            a = float(a)
        except ValueError:
            raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.single_slit_diffraction_solve:\n"
            "\t\tThird input must be a floating point number.\n"
            "\t\tYour input has type: %s\n"
            % (type(a).__name__)
            )

    f = np.sinc(a*x/z)*np.sinc(a*x/z)
    return f

def double_slit_diffraction(x, z, a, d):
    """
        Purpose:
            Compute diffraction through a double slit for the
            variable x with a distance z from the slit and
            slit parameter a and a distance d between the slits.
            This assumes Fraunhofer diffraction.
        Variables:
            :x:
                A real or complex argument, or numpy array.
            :z (*float*):
                The perpendicular distance from the slit
                plane to the observer.
            :a (*float*):
                The slit parameter. This is a unitless paramter
                defined as the ratio between the slit width and
                the wavelength of the incoming signal.
            :d (*float*):
                The distance between slits.
        Outputs:
            :f:
                Single slit diffraction pattern.
        Dependences:
            [1] numpy
    """
    y = x
    try:
        x = np.array(x)
        if (not np.all(np.isreal(x))) and (not np.all(np.iscomplex(x))):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\tspecial_function.double_slit_diffraction_solve:\n"
                "\t\tInput must be a real or complex valued numpy array.\n"
                "\t\tThe elements of your array have type: %s"
                % (x.dtype)
            )
        else:
            del y
    except (TypeError, ValueError) as errmes:
        raise TypeError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.double_slit_diffraction_solve:\n"
            "\t\tFirst input must be a real or complex valued numpy array.\n"
            "\t\tYour input has type: %s\n"
            "\tOriginal Error Mesage: %s\n"
            % (type(y).__name__, errmes)
        )
    if not (isinstance(z, float)):
        try:
            z = float(z)
        except ValueError:
            raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.double_slit_diffraction_solve:\n"
            "\t\tSecond input must be a floating point number.\n"
            "\t\tYour input has type: %s\n"
            % (type(z).__name__)
            )
    if not (isinstance(a, float)):
        try:
            a = float(a)
        except ValueError:
            raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.double_slit_diffraction_solve:\n"
            "\t\tThird input must be a floating point number.\n"
            "\t\tYour input has type: %s\n"
            % (type(a).__name__)
            )
    if not (isinstance(d, float)):
        try:
            d = float(d)
        except ValueError:
            raise ValueError(
            "\n\tError Encountered:\n"
            "\trss_ringoccs: Diffcorr Subpackage\n"
            "\tspecial_function.double_slit_diffraction_solve:\n"
            "\t\tFourth input must be a floating point number.\n"
            "\t\tYour input has type: %s\n"
            % (type(a).__name__)
            )
    f1 = np.sinc(a*x/z)*np.sinc(a*x/z)
    f2 = np.sin(TWO_PI*d*x/z)*np.sin(TWO_PI*d*x/z)
    f3 = 4.0*np.sin(ONE_PI*d*x/z)*np.sin(ONE_PI*d*x/z)
    f = f1*f2/f3

    return f

def sq_well_solve(x, a, b, F, invert=False):
    """
        Function:
            sq_well_solve
        Purpose:
            Computes the solution of diffraction through a square well.
        Variables:
            :x:
                Real numpy array. The independent variable.
            :a (*float*):
                The LEFTMOST endpoint of the square well.
            :b (*float*):
                The RIGHTMOST endpoint of the square well.
            :F (*float*):
                The Fresnel scale.
        Output:
            :H:
                Complex numpy array.
                Diffraction pattern of a square well on
                the interal [a,b].
        History:
            Translated from IDL: RJM - 2018/05/15 8:03 P.M.
            Updated and added error checks: RJM - 2018/09/19 7:19 P.M.
    """
    if (not isinstance(x, np.ndarray)):
        raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec subpackage\n"
                "\tspecial_functions: sq_well_solve\n"
                "\tFirst variable should be a numpy array.\n"
        )
    elif (not np.isreal(x).all()):
        raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec subpackage\n"
                "\tspecial_functions: sq_well_solve\n"
                "\tFirst variable should be real valued.\n"
        )
    else:
        pass

    if (not isinstance(a,float)):
        try:
            a = float(a)
        except TypeError:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec subpackage\n"
                "\tspecial_functions: sq_well_solve\n"
                "\tSecond variable should be a floating point number.\n"
                "\tYour input has type: %s"
                % (type(a).__name__)
            )
    else:
        pass

    if (not isinstance(b,float)):
        try:
            b = float(b)
        except TypeError:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec subpackage\n"
                "\tspecial_functions: sq_well_solve\n"
                "\tThird variable should be a floating point number.\n"
                "\tYour input has type: %s"
                % (type(b).__name__)
            )
    else:
        pass

    if (not isinstance(F,float)):
        try:
            F = float(F)
        except TypeError:
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec subpackage\n"
                "\tspecial_functions: sq_well_solve\n"
                "\tFourth variable should be a floating point number.\n"
                "\tYour input has type: %s"
                % (type(b).__name__)
            )
    else:
        pass

    H = (0.5 - 0.5j) * (fresnel_cos((b-x)/F)-fresnel_cos((a-x)/F)+
                        1j*(fresnel_sin((b-x) / F)-fresnel_sin((a-x)/F)))

    if not invert:
        H = 1-H

    return H
