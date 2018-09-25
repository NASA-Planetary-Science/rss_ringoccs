import numpy as np
from scipy.special import erf, lambertw
from . import window_functions as wf
import pdb
SQRT_PI_2 = 0.886226925452758013649084

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
        Function:
            savitsky_golay
        Purpose:
            To smooth data with a Savitzky-Golay filter.
            This removes high frequency noise while
            maintaining many of the original features of
            the input data.
        Parameters
        ----------
        y : array_like, shape (N,)
            the values of the time history of the signal.
        window_size : int
            the length of the window. Must be an odd integer number.
        order : int
            the order of the polynomial used in the filtering.
            Must be less then `window_size` - 1.
        deriv: int
            the order of the derivative to compute (default = 0 means only smoothing)
        Returns
        -------
        ys : ndarray, shape (N)
            the smoothed signal (or it's n-th derivative).
        Notes
        -----
        The Savitzky-Golay is a type of low-pass filter, particularly
        suited for smoothing noisy data. The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        Examples
        --------
        t = np.linspace(-4, 4, 500)
        y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
        ysg = savitzky_golay(y, window_size=31, order=4)
        import matplotlib.pyplot as plt
        plt.plot(t, y, label='Noisy signal')
        plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
        plt.plot(t, ysg, 'r', label='Filtered signal')
        plt.legend()
        plt.show()
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
        Data by Simplified Least Squares Procedures. Analytical
        Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Cambridge University Press ISBN-13: 9780521880688
    """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")

    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    elif window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    else:
        pass

    half_window = (window_size -1) // 2

    # precompute coefficients
    b = np.zeros((window_size, order+1))
    b[..., 0] = 1
    for k in range((window_size - 1) // 2):
        n0 = ((window_size-1) // 2) - k
        m = n0
        n = -n0
        for j in range(1,order+1):
            b[k, j] = n
            b[window_size-1-k, j] = m
            n *= -n0
            m *= n0

    b = np.mat(b)

    m = 1
    for i in range(deriv):
        m *= rate*(deriv-i)

    m *= np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

def fresnel_transform(rho_km_vals, F_km_vals, phi_rad_vals, B_rad_vals,
                    D_km_vals, T_vals, lambda_km_vals, w_km_vals,
                    dx_km, wtype, start, n_used, norm=True, fft=False,
                    verbose=True, psitype='full'):
    """
        Procedure: fresnel_forward
        Purpose:   Computes the forward model of diffraction from a set of
            reconstructed data.
        Variables:
            rho_vals:        Ring radius, in kilometers.
            F_vals:          Fresnel scale, in kilometers.
            phi_rad_vals:    Ring azimuth angle, in radians.
            B_rad_vals:      Ring opening angle, in radians.
            lambda_sky_vals: Wavelength of recieved signal, in kilometers.
            D_vals:          Spacecraft-RIP distance, in kilometers.
            dx:              Sampling spacing, in kilometers.
            T_vals:          Reconstructed complex transmittance.
            w_vals:          Window width, in kilometers.
            wtype:           Window used in reconstruction, string.
            start:           Starting point of reconstructed data.
            n_used:          Number of reconstructed points.
        Keywords:
            Normalize: Parameter for normalizing the complex transmittance by
                the window function that is used. Default is True. Set to False
                to skip this feature.
        Output:
            phase_fwd_vals  - Phase of the forward model, in radians.
            T_hat_fwd_vals  - Complex transmittance of forward model.
            p_norm_fwd_vals - Normalized power of forward model, unitless.
        History:
            Translated from IDL: RJM - 2018/05/14 5:06 P.M.
    """
    w_max     = np.max(w_km_vals[start:start + n_used])
    nw_fwd    = int(np.ceil(w_max / (2.0 * dx_km)))
    start_fwd = int(start + nw_fwd)
    n_fwd     = int(n_used - 2 * nw_fwd)

    # Compute product of wavenumber and RIP distance.
    kD_vals = 2.0 * np.pi * D_km_vals / lambda_km_vals

    # Compute Cosine of opening angle.
    cosb = np.cos(B_rad_vals)

    # Precompute variables for speed.
    cosbD = cosb/D_km_vals
    cosb2 = cosb*cosb

    # Compute cosine and sine of ring azimuth angle.
    cosphi0 = np.cos(phi_rad_vals)
    sinphi0 = np.sin(phi_rad_vals)

    # Precompute squares of variables for speed.
    dsq = D_km_vals*D_km_vals
    rsq = rho_km_vals*rho_km_vals

    # Define window function.
    fw = wf.func_dict[wtype]["func"]

    # Define normalization function.
    nrm = wf.normalize

    # Set inverse function to FFT or Integration.
    if fft:
        finv = fresnel_inverse_fft
    else:
        finv = fresnel_inverse
    # Set psi approximation.
    if (psitype == 'full'):
        psif = psi_func
    elif (psitype == 'mtr4'):
        psif = psi_quartic
    else:
        psif = psi_func

    # Create an empty array for reconstructed complex transmittance.
    T_hat_fwd_vals = 0.0 * T_vals

    # Set the first computed point to the 'start' variable.
    center = start_fwd

    # Compute first window width and window function.
    w_init = w_km_vals[center]
    w_func = fw(w_init, dx_km, error_check=False)

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
        for i in np.arange(n_fwd):

            # Current point being computed.
            center = start_fwd+i

            # Ring radius of current point.
            r0 = rho_km_vals[center]

            # Window width for current point.
            w = w_km_vals[center]

            # Window function for current point.
            w_func = fw(w, dx_km, error_check=False)

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
            T = T_vals[crange]

            # Compute 'approximate' Fresnel Inversion for current point
            T_hat_fwd_vals[center] = finv(T_hat, ker, dx_km, F, error_check=False)

            # If normalization has been set, normalize the reconstruction
            if norm:
                T_hat_fwd_vals[center] *= nrm(dx_km, ker, F)
            if verbose:
                print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d Inversion"
                      % (i, n_fwd, nw, loop), end="\r")
    else:
        for i in np.arange(n_fwd):

            # Current point being computed.
            center = start_fwd+i

            # Compute current radius and RIP distance.
            r0 = rho_km_vals[center]
            d2 = dsq[center]

            # Compute square of ring radius.
            r02 = rsq[center]
            phi0 = phi_rad_vals[center]

            # Precomputed variables for speed.
            cbd = cosbD[center]
            cb2 = cosb2[center]
            cp0 = cosphi0[center]
            sp0 = sinphi0[center]

            # Compute product of wavenumber and RIP Distance.
            kD = kD_vals[center]

            # Current window width and Fresnel scale.
            d = D_km_vals[center]
            b = B_rad_vals[center]
            w = w_km_vals[center]
            F = F_km_vals[center]

            # If the window width has changed, recompute variables.
            if (np.abs(w_init - w) >= 2.0 * dx_km):

                # Reset w_init and recompute window function.
                w_init = w
                w_func = fw(w, dx_km, error_check=False)

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
            psi_vals = kD * psif(r, r0, phi_s_rad, phi0, d, b, error_check=False)

            # Compute kernel function for Fresnel inverse
            ker = w_func*np.exp(1j*psi_vals)

            # Range of diffracted data that falls inside the window
            T = T_vals[crange]

            # Compute 'approximate' Fresnel Inversion for current point
            T_hat_fwd_vals[center] = finv(T, ker, dx_km, F, error_check=False)

            # If normalization has been set, normalize the reconstruction
            if norm:
                T_hat_fwd_vals[center] *= nrm(dx_km, ker, F)
            if verbose:
                print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d Forward"
                      % (i, n_fwd, nw, loop), end="\r")
    return T_hat_fwd_vals

def compute_norm_eq(w_func):
    """
        Function:
            compute_norm_eq
        Purpose:
            Compute normalized equivalenth width of a given function.
        Variables:
            w_func: Any function (Usually a window function).
        Outputs:
            normeq: The normalized equivalent width of w_func.
        Dependencies:
            [1] sys
            [2] numpy
        Notes:
            The normalized equivalent width is effectively computed
            using Riemann sums to approximate integrals. Therefore
            large dx values (Spacing between points in w_func)
            will result in an inaccurate normeq. One should keep
            this in mind during calculations.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        Examples:
            Use diffcorr to compute the Kaiser-Bessel 2.5 window
            of witdh 30 and spacing 0.1, and then use compute_norm_eq
            to compute the normalized equivalent width:
                In [1]: import diffcorr as dc
                In [2]: w = dc.kb25(30,0.1)
                In [3]: normeq = dc.compute_norm_eq(w)
                In [4]: print(normeq)
                1.6573619266424229
            In contrast, the actual value is 1.6519208.
            Compute the normalized equivalent width for the squared
            cosine window of width 10 and spacing 0.25.
                In [1]: import diffcorr as dc
                In [2]: w = dc.coss(10,0.25)
                In [3]: normeq = dc.compute_norm_eq(w)
                In [4]: print(normeq)
                1.5375000000000003
            The normalized equivalent width of the squared cosine
            function can be computed exactly using standard methods
            from a calculus course. It's value is exactly 1.5
            If we use a smaller dx when computing w, we get a better
            approximation. Use width 10 and spacing 0.001.
                In [1]: import diffcorr as dc
                In [2]: w = dc.coss(10,0.001)
                In [3]: normeq = dc.compute_norm_eq(w)
                In [4]: print(normeq)
                1.50015
        History:
            Created: RJM - 2018/05/16 3:54 P.M.
    """
    if not check_real(w_func):
        raise TypeError("Input must be real valued")
    nw      = np.size(w_func)
    normeq  = nw*(np.sum(w_func**2)) / (np.sum(w_func)**2)
    return normeq

def fresnel_scale(Lambda,d,phi,b,DEG=False):
    """
    Function:
        fresnel_scale
    Purpose:
        Compute the Fresnel Scale from lambda, D, Phi, and B.
    Variables:
        Lambda: Wavelength of the incoming signal.
        d:      RIP-Spacecraft Distance.
        phi:    Ring Azimuth Angle.
        b:      Ring Opening Angle.
    Keywords:
        DEG:    Set True if phi/b are in degrees. Default is radians.
    Output:
        FRES:   The Fresnel scale.
    NOTE:
        Lambda and d must be in the same units. The output (Fresnel
        scale) will have the same units as lambda and d. In addition,
        b and phi must also have the same units. If b and phi are in
        degrees, make sure to set DEG = True. Default is radians.
    History:
        Translated from IDL: RJM - 2018/04/15 12:36 P.M.
    """
    if (not check_real(Lambda)):
        sys.exit("Lambda must be real")
    if (not check_real(d)):
        sys.exit("D must be real")
    if (not check_real(phi)):
        sys.exit("Phi must be real")
    if (not check_real(b)): 
        sys.exit("B must be real")
    if DEG:
        cb = np.cos(b * np.pi / 180)
        sb = np.sin(b * np.pi / 180)
        sp = np.sin(phi * np.pi / 180)
    else:
        cb = np.cos(b)
        sb = np.sin(b)
        sp = np.sin(phi)
    fres = np.sqrt(0.5 * Lambda * d * (1 - (cb**2) * (sp**2)) / (sb**2))
    return fres

def psi_d1_phi(r,r0,d,b,phi,phi0,error_check=True):
    """
        Function: psi_d1_phi
        Purpose:  Calculate dpsi/dphi from geometry variables.
        Variables:
            r:    Ring radius variable (Integrated over). km.
            r0:   Ring intercept point. km.
            d:    RIP-Spacecraft distance. km.
            b:    The ring opening angle.
            phi:  The ring azimuth angle (variable).
            phi0: The ring azimuth angle (value).
        History:
            Translated from IDL: RJM - 2018/05/15 7:06 P.M.
    """
    if error_check:
        if (not check_real(r)):
            sys.exit("r must be real valued")
        if (not check_real(r0)):
            sys.exit("r0 must be real valued")
        if (not check_real(d)):
            sys.exit("d must be real valued")
        if (not check_real(b)):
            sys.exit("b must be real valued")
        if (not check_real(phi)):
            sys.exit("phi must be real valued")
        if (not check_real(phi0)):
            sys.exit("phi0 must be real valued")
    else: pass
    cb   = np.cos(b)
    sp   = np.sin(phi)
    cp   = np.cos(phi)
    sp0  = np.sin(phi0)
    cp0  = np.cos(phi0)
    xi   = (cb / d) * (r0*cp0 - r*cp)
    eta  = (r0**2 + r**2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / (d**2)
    v1   = r * cb * sp / d
    v2   = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d**2)
    psi_d1_phi_vals = (2.0*v1 + v2) / (2.0 * np.sqrt(1.0 + 2.0*xi + eta)) - v1
    return psi_d1_phi_vals

def psi_d2_phi(r,r0,d,b,phi,phi0,error_check=True):
    """
        Function: psi_d2_phi_fast
        Purpose:  Calculate second derivative of psi with respect to phi from
        geometry variables using previously computed sines and cosines.
        Variables:
            r:    Ring radius variable (Integrated over). km.
            r0:   Ring intercept point. km.
            d:    RIP-Spacecraft distance. km.
            cb:   The cosine of the ring opening angle.
            cp:   The cosine of the ring azimuth angle (variable).
            sp:   The sine of the ring azimuth angle (variable).
            cp0:  The cosine of the ring azimuth angle (value).
            sp0:  The sine of the ring azimuth angle (value).
        History:
            Translated from IDL: RJM - 2018/05/15 5:36 P.M.
    """
    if error_check:
        if (not check_real(r)):
            sys.exit("r must be real valued")
        if (not check_real(r0)):
            sys.exit("r0 must be real valued")
        if (not check_real(d)):
            sys.exit("d must be real valued")
        if (not check_real(b)):
            sys.exit("b must be real valued")
        if (not check_real(phi)):
            sys.exit("phi must be real valued")
        if (not check_real(phi0)):
            sys.exit("phi0 must be real valued")
    else: pass
    cb   = np.cos(b)
    sp   = np.sin(phi)
    cp   = np.cos(phi)
    sp0  = np.sin(phi0)
    cp0  = np.cos(phi0)
    xi    = (cb / d) * (r0*cp0 - r*cp)
    eta   = ((r0**2) + (r**2) - 2.0*r*r0*(sp*sp0 + cp*cp0)) / (d**2)
    v1    = r * cb * cp / d
    v2    = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / (d**2)
    v3    = r * cb * sp / d
    v4    = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d**2)
    dphia = (2.0*v1 + v2)/(2.0 * np.sqrt(1.0 + 2.0*xi + eta))
    dphib = v1 + ((2.0*v3 + v4)**2)/(4.0 * (np.sqrt(1.0 + 2.0*xi + eta)**3))
    psi_d2_phi_vals = dphia - dphib
    return psi_d2_phi_vals

def fresnel_inverse(T,ker,DX,f_scale):
    """
        Function: fresnel_transform
        Purpose:  Compute the approximate inverse of a Fresnel transform.
            Variables:
                T:   The input function (Diffraction data).
                KER: The Fresnel kernel.
                dx:  The size of a bin (r[1]-r[0]).
                f_scale: The Frensel scale. If not set, a default of 1 is used.
        Output:
            T_HAT:   The forward model for the diffraction pattern.
        History:
            Translated from IDL: RJM - 2018/04/15 12:23 P.M.
    """
    if (not check_real(T)) and (not check_complex(T)):
        sys.exit('T must be real or complex')
    if (not check_real(ker)) and (not check_complex(ker)):
        sys.exit('Ker must be real or complex')
    if (np.size(T) != np.size(ker)):
        sys.exit('T and ker have a different number of elements')
    if (not check_pos_real(DX)):
        sys.exit('DX must be a positive number')
    if (not check_pos_real(f_scale)):
        sys.exit('f_scale must be positive')
    T_hat = np.sum(ker * T) * DX * (1.0-1.0j) / (2. * f_scale)
    return T_hat

def fresnel_inverse_fft(T_hat,ker,dx,f_scale,error_check=True):
    if error_check:
        if (not check_real(T_hat)) and (not check_complex(T_hat)):
            sys.exit('T_hat must be real or complex')
        if (not check_real(ker)) and (not check_complex(ker)):
            sys.exit('Ker must be real or complex')
        if (np.size(T_hat) != np.size(ker)):
            sys.exit('T_hat and ker have a different number of elements')
        if (not check_pos_real(dx)):
            sys.exit('DX must be a positive number')
        if (not check_pos_real(f_scale)):
            sys.exit('F_SCALE must be positive')
    else: pass
    nw = np.size(T_hat)
    fft_t_hat       = np.fft.fft(T_hat)
    fft_conv        = np.fft.fft(ker)
    inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
    inv_t_hat      *= dx*(np.complex(1.0,1.0))/(2.0*f_scale)
    T               = inv_t_hat[int((nw-1)/2)]
    return T

def psi_func(r,r0,d,b,phi,phi0,error_check=True):
    """
        Function: psi_fast
        Purpose:  Calculate psi from geometry variables.
        Variables:
            r:    Ring radius variable.
            r0:   Ring intercept point.
            D:    RIP-Spacecraft distance.
            b:    Ring opening angle.
            phi:  Ring azimuth angle (variable).
            phi0: Ring azimuth angle (variable).
        History:
            Translated from IDL: RJM - 2018/05/15 7:48 P.M.
    """
    if error_check:
        if (not check_real(r)):
            sys.exit("r must be real valued")
        if (not check_real(r0)):
            sys.exit("r0 must be real valued")
        if (not check_real(d)):
            sys.exit("d must be real valued")
        if (not check_real(b)):
            sys.exit("b must be real valued")
        if (not check_real(phi)):
            sys.exit("phi must be real valued")
        if (not check_real(phi0)):
            sys.exit("phi0 must be real valued")
    else: pass
    cb   = np.cos(b)
    sp   = np.sin(phi)
    cp   = np.cos(phi)
    sp0  = np.sin(phi0)
    cp0  = np.cos(phi0)
    xi   = (cb / d) * (r0*cp0 - r*cp)
    eta  = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d*d)
    psi_vals   = np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi)
    return psi_vals

def resolution_inverse(x):
    """
        Function:
            resolution_inverse
        Purpose:
            Compute the inverse of y = x/(exp(-x)+x-1)
        Variables:
            x:
                A real or complex number, or numpy array.
        Outputs:
            f:
                The inverse of x/(exp(-x)+x-1)
        Dependencies:
            [1] numpy
            [2] scipy.special
        Method:
            The inverse of x/(exp(-x)+x-1) is computed using the
            LambertW function. This function is the inverse of
            y = x * exp(x). This is computed using the scipy.special
            subpackage using their lambertw function.
        Warnings:
            [1] The real part of the argument must be greater than 1.
            [2] The scipy.special lambertw function is slightly
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
        History:
            Translated from IDL: RJM - 2018/05/15 1:49 P.M.
            Added error checks: RJM - 2018/09/19 6:57 P.M.
    """
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
        P1 = x/(1.0-x)
        P2 = P1*np.exp(P1)
        f = lambertw(P2)-P1

    if np.all(np.isreal(x)):
        f = np.real(f)

    return f

def fresnel_cos(x):
    """
        Function:
            fresnel_cos
        Purpose:
            Compute the Fresnel cosine function.
        Variables:
            x:
                A real or complex number, or numpy array.
        Outputs:
            f_cos:
                The fresnel cosine integral of x.
        Dependences:
            [1] numpy
            [2] scipy
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
    f_cos = ((0.25-0.25j)*erf((1.0+1.0j)*x*SQRT_PI_2)+
            (0.25+0.25j)*erf((1.0-1.0j)*x*SQRT_PI_2))
    if (np.isreal(x).all()):
        f_cos = np.real(f_cos)
    return f_cos

def fresnel_sin(x):
    """
        Function:
            fresnel_sin
        Purpose:
            Compute the Fresnel sine function.
        Variables:
            x:
                A real or complex argument, or numpy array.
        Outputs:
            f_sin:
                The fresnel sine integral of x.
        Dependences:
            [1] numpy
            [2] scipy
        Notes:
            [1] The Fresnel sine integral is the solution to the
                equation dy/dx = sin(pi/2 * x^2), y(0) = 0. In other
                words, y = integral (t=0 to x) sin(pi/2 * t^2) dt
            [2] The Fresnel Cossine and Sine integrals are computed by
                using the scipy.special Error Function. The Error
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
        History:
            Translated from IDL: RJM - 2018/05/14 3:53 P.M.
            Allow complex arguments: RJM - 2018/05/16 1:26 P.M.
            Updated and added error check: RJM - 2018/09/19 7:01 P.M.
    """
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
    f_sin = ((0.25+0.25j)*erf((1.0+1.0j)*x*SQRT_PI_2)+
            (0.25-0.25j)*erf((1.0-1.0j)*x*SQRT_PI_2))
    if (np.isreal(x).all()):
        f_sin = np.real(f_sin)
    return f_sin

def single_slit_diffraction_solve(x, z, a):
    """
        Function:
            single_slit_diffraction_solve
        Purpose:
            Compute diffraction through a single slit for the
            variable x with a distance z from the slit and
            slit parameter a. This assume Fraunhofer diffraction.
        Variables:
            x:
                A real or complex argument, or numpy array.
            z:
                Float
                The perpendicular distance from the slit plane to
                the observer.
            a:
                Float
                The slit parameter. This is a unitless paramter
                defined as the ratio between the slit width and
                the wavelength of the incoming signal.

        Outputs:
            f:
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

def double_slit_diffraction_solve(x, z, a, d):
    """
        Function:
            double_slit_diffraction_solve
        Purpose:
            Compute diffraction through a double slit for the
            variable x with a distance z from the slit and
            slit parameter a and a distance d between the slits.
            This assume Fraunhofer diffraction.
        Variables:
            x:
                A real or complex argument, or numpy array.
            z:
                Float
                The perpendicular distance from the slit plane to
                the observer.
            a:
                Float
                The slit parameter. This is a unitless paramter
                defined as the ratio between the slit width and
                the wavelength of the incoming signal.
            d:
                Float
                The distance between slits.

        Outputs:
            f:
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
    f2 = np.sin(2.0*np.pi*d*x/z)*np.sin(2.0*np.pi*d*x/z)
    f3 = 4.0*np.sin(np.pi*d*x/z)*np.sin(np.pi*d*x/z)
    f = f1*f2/f3
    return f

def sq_well_solve(x, a, b, F, invert=False):
    """
        Function:
            sq_well_solve
        Purpose:
            Computes the solution of diffraction through a square well.
        Variables:
            x:
                Real numpy array.
                The independent variable.
            a:
                Float
                The LEFTMOST endpoint of the square well.
            b:
                Float
                The RIGHTMOST endpoint of the square well.
            F: The Fresnel scale.
        Output:
            H:
                Complex numpy array.
                Diffraction pattern of a square well on
                the interal [a,b].
        History:
            Translated from IDL: RJM - 2018/05/15 8:03 P.M.
            Updated and added error checks: RJM - 2018/09/19 7:19 P.M.
    """
    if (np.size(a) != 1):
        raise TypeError("Endpoints must be floating point numbers")
    elif (not isinstance(a,float)):
        try:
            a = float(a)
        except TypeError:
            raise TypeError("Endpoints must be floating point numbers")
    elif (not np.isreal(a)):
        raise ValueError("Endpoints must be real valued")
    if (np.size(b) != 1):
        raise TypeError("Endpoints must be floating point numbers")
    elif (not isinstance(b,float)):
        try:
            b = float(b)
        except TypeError:
            raise TypeError("Endpoints must be floating point numbers")
    elif (not np.isreal(b)):
        raise ValueError("Endpoints must be real valued")
    if (np.size(F) != 1):
        raise TypeError("Endpoints must be floating point numbers")
    elif (not isinstance(F,float)):
        try:
            F = float(F)
        except TypeError:
            raise TypeError("Endpoints must be floating point numbers")
    elif (not np.isreal(F)):
        raise ValueError("Endpoints must be real valued")
    if (not isinstance(x, np.ndarray)):
        raise TypeError("Independant variable must be a numpy array")
    elif (not np.isreal(x).all()):
        raise ValueError("Independent variable must be real valued")
    H = (0.5 - 0.5j) * (
        fresnel_cos((b - x) / F) - fresnel_cos((a - x) / F)
        + 1j*(fresnel_sin((b - x) / F) - fresnel_sin((a - x) / F))
    )
    if not invert:
        H = 1-H
    return H
