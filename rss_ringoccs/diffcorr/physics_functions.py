import numpy as np
from .window_functions import rect, coss, kb20, kb25, kb35, kbmd20, kbmd25
from .window_functions import normalize, func_dict
from scipy.constants import speed_of_light

def power_func(T_in):
    """
        Function:
            power
        Purpose:
            Compute power from complex transmittance.
        Variables:
            T_in:   The complex transmittance.
        Output:
            power:  The power.
        Dependencies:
            [1] numpy
            [2] sys
        Notes:
            [1] The power is simply the square of the absolute value
                of the complex transmittance. This equation works for
                both a diffracted and a reconstructed transmittance.
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
        Examples:
            Compute the power of a diffraction pattern through a
            square well.
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01-50.0
                In [5]: T = dc.sq_well_solve(x,-5,5,1)
                In [6]: p = dc.power_func(T)
                In [7]: plt.show(plt.plot(x,p))
        History:
            Created: RJM - 2018/05/16 5:19 A.M.
    """
    if (not check_real(T_in)) and (not check_complex(T_in)):
        raise TypeError("Complex transmittance must be real or complex valued.")
    else:
        power = (np.abs(T_in))**2
    return power

def phase_func(T_in):
    """
        Function:
            phase_func
        Purpose:
            Compute the phase from the complex transmittance.
        Variables:
            T_in:T  The complex transmittance.
        Output:
            phase:  The phase (in radians).
        Dependencies:
            [1] numpy
            [2] sys
        Notes:
            [1] The phase of the complex transmittance is the angle
                made it makes with the x-axis in the complex plane.
                This is the arctangent of the ratio of the imaginary
                part to the real part. This equation works for both
                diffracted and reconstructed transmittances.
        Refernces:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] S. W. Asmar, R. G. French, E. A. Marouf, P. Schinder,
                J. W. Armstrong, P. Tortora, L. Iess, A. Anabtawi,
                A. J. Kliore, M. Parisi, M. Zannoni, and D. Kahan,
                Cassini Radio Science User's Guide, September, 2014,
                https://pds-rings.seti.org/cassini/rss/
        Examples:
            Calculate the phase from the diffraction through a well.
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01-50.0
                In [5]: T = dc.sq_well_solve(x,-5,5,1)
                In [6]: phi = dc.phase_func(T)
                In [7]: plt.show(plt.plot(x,phi))
        History:
            Created: RJM - 2018/05/16 5:19 A.M.
    """
    if (not check_real(T_in)) and (not check_complex(T_in)):
        raise TypeError("Complex transmittance must be real or complex valued.")
    else:
        phase = np.arctan2(np.imag(T_in),np.real(T_in))
    return phase

def tau_func(T_in,mu):
    """
        Function:
            tau_func
        Purpose:
            Compute the normalized optical depth.
        Variables:
            T_in:   The complex transmittance.
            mu:     The sine of the ring opening angle.
        Output:
            tau:    The normalized optical depth.
        Dependencies:
            [1] numpy
            [2] sys
        Notes:
            [1] The optical depth is proportional to the natural
                logarithm of the transmitted power. It is a unitless
                variable that arises in the study of radiative
                transfer, in particular the transfer equation. The
                normalized optical depth is normalized version of
                this, taking geometrical factors into account and
                using the normalized power.
            [2] Tau is often used to represent the optical depth. The
                equation used it tau = mu * ln(Power), where mu is
                the sine of the ring opening angle (Denoted B), and
                where ln is the natural logarithm.
        References:
            [1] George B. Rybicki and Alan P. Lightman,
                Radiative Processes in Astrophysics,
                Wiley, 29 December 2007
            [2] S. Chandrasekhar, Radiative Transfer,
                Dover Publications, 1960
            [3] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [4] S. W. Asmar, R. G. French, E. A. Marouf, P. Schinder,
                J. W. Armstrong, P. Tortora, L. Iess, A. Anabtawi,
                A. J. Kliore, M. Parisi, M. Zannoni, and D. Kahan,
                Cassini Radio Science User's Guide, September, 2014,
                https://pds-rings.seti.org/cassini/rss/
        Example:
            Plot the normalized optical depth of transmittance
            through a square well using a ring opening angle of
                B = pi/2 (Or, mu = 1).
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: T = dc.sq_well_solve(x,-5,5,1)
                In [7]: tau = dc.tau_func(T,1)
                In [8]: plt.show(plt.plot(x,tau))
        History:
            Created: RJM - 2018/05/16 7:18 A.M.
    """
    if (not check_real(T_in)) and (not check_complex(T_in)):
        raise TypeError("Complex transmittance must be real or complex valued.")
    elif (not check_real(mu)):
        raise TypeError("mu must be real valued.")
    else:
        p           = power_func(T_in)
        crange      = (p>0).nonzero()
        tau         = np.zeros(np.size(p))
        tau[crange] = -mu[crange] * np.log(np.abs(p[crange]))
    return tau

def wker(w,psi):
    """
        Function:
            wker
        Purpose:
            Compute the weighted kernel function.
        Variables:
            psi:  The independent variable.
            w:    The weight function.
        Output:
            kernel: The weighted kernel function.
        History:
            Translated from IDL: RJM - 2018/05/16 5:10 A.M.
    """
    kernel = w * np.exp(1j * psi)
    return kernel

def freq_wav(freqwav):
    """
        Function:
            freq_wav
        Purpose:
            Converts frequency to wavelength, and vice versa.
        Variables:
            FREQWAV: Frequency (wavelength) of the input in Hz (km).
        Outputs:
            WAVFREQ: Wavelength (frequency) of the input in Km (Hz).
        NOTE:
            Frequency MUST be Hertz, wavelength MUST be Kilometers.
        History:
            Translated from IDL: RJM - 2018/05/14 11:41 A.M.
    """
    if not check_real(freqwav):
        sys.exit("Input must be real valued")
    elif (np.min(freqwav) <= 0):
        sys.exit("Input must be positive")
    else:
        wavfreq = speed_of_light*0.001 / freqwav
    return wavfreq

def fresnel_forward(rho_km_vals,F_km_vals,phi_rad_vals,B_rad_vals,D_km_vals,
    T_vals,lambda_km_vals,w_km_vals,dx_km,wtype,start,n_used,norm=True,fft=False,
    verbose=True,psitype='full'):
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
    fw = func_dict[wtype]["func"]

    # Define normalization function.
    nrm = normalize

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

def fresnel_inversion(rho_vals,F_vals,phi_rad_vals,B_vals,D_vals,
    T_hat_vals,lambda_vals,w_vals,dx,wtype,start,n_used,norm=True,fft=False,
    verbose=True,psitype='full'):
    """
        Function:  fresnel_inversion
        Purpose:   Computes the fresnel inversion from a set of diffracted data
            using a 'fast' method to speed up computation time. This is
            achieved by computing cosine and sine function in the outer for
            loop, and then passing these computed values into the functions
            that need them. The normal version passes the arguments to the
            functions, and then cosines and sines are computed within the
            function. For small data sets or coarser resolutions, the normal
            version is faster. Both the normal and fast versions output
            completely identical results.
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
            T_vals  - Reconstructed Complex transmittance.
        History:
            Translated from IDL: RJM - 2018/05/16 6:26 A.M.
    """
    # Compute necessary variables.
    kD_vals   = 2.0*np.pi*D_vals/lambda_vals
    # Define functions
    fw        = func_dict[wtype]["func"]
    nrm       = normalize
    if fft:
       finv   = fresnel_inverse_fft
    else: 
        finv  = fresnel_inverse
    # Calculate the corrected complex amplitude, point by point
    T_vals    = T_hat_vals * 0.0
    w_init    = w_vals[start]
    w_func    = fw(w_init,dx)
    nw        = np.size(w_func)
    crange    = np.arange(int(start-(nw-1)/2),int(1+start+(nw-1)/2))
    phi_s_rad1 = phi_rad_vals[start]
    if psitype == 'taylor2':
        for i in np.arange(n_used):
            center      = start+i
            r0          = rho_vals[center]
            w           = w_vals[center]
            w_func      = fw(w,dx,error_check=False)
            nw          = np.size(w_func)
            crange      = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r           = rho_vals[crange]
            F           = F_vals[center]
            x           = (r-r0)/F
            psi_vals    = (np.pi/2.0)*x*x
            ker         = w_func*np.exp(-1j*psi_vals)
            T_hat       = T_hat_vals[crange]
            
            T_vals[center] = finv(T_hat,ker,dx,F)
            if norm:T_vals[center] *= nrm(r,w_func,F,error_check=False)
            if verbose:
                print("Pt: %d  Tot: %d  Width: %d  Inversion"\
                % (i,n_used,nw),end="\r")
    elif (psitype == 'mtr2') or (psitype == 'mtr3') or (psitype == 'mtr4'):
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            d      = D_vals[center]
            b      = B_vals[center]
            phi0   = phi_rad_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            F      = F_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx,error_check=False)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                dphi_s_rad = psi_factor(r,r0,b,phi0,error_check=False)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
            else:
                crange    += 1
                r          = rho_vals[crange]
                phi_s_rad  = phi_s_rad1
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
            loop = 0
            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                loop      += 1
                if loop > 5:
                    break
            phi_s_rad1 = phi_s_rad

            # Compute psi and then compute averages and differences across psi.
            psi_full   = kD * psi(r,r0,d,b,phi_s_rad,phi0,error_check=True)
            n1         = 0                                 #Left Endpoint
            n2         = np.min((r0-w/4<=r).nonzero())     #Left midpoint
            n3         = np.max((r0+w/4>=r).nonzero())     #Right midpoint
            n4         = nw-1                              #Right endpoint
            d_psi_half = psi_full[n3]-psi_full[n2]         #Midpoint difference
            d_psi_full = psi_full[n4] - psi_full[n1]       #Endpoint difference
            a_psi_half = (psi_full[n3]+psi_full[n2])/2     #Midpoint average
            a_psi_full = (psi_full[n1]+psi_full[n4])/2     #Endpoint average
            x          = (r-r0)
            w          = np.max(x)-np.min(x)

            #Compute coefficients for the polynomial expansion
            c1  = (8.0*d_psi_half-d_psi_full)/(3.0*w)           #Linear term
            c2  = 4.0*(16.0*a_psi_half-a_psi_full)/(3.0*w*w)    #Quadratic term
            c3  = 16.0*(d_psi_full-2.0*d_psi_half)/(3.0*w*w*w)  #Cubic term
            c4  = 64.0*(a_psi_full-4.0*a_psi_half)/(3.0*w*w*w*w)#Quartic term

            psi_vals = c1*x + c2*x*x        #Second order appoximation
            if (psitype == 'mtr3'):
                psi_vals += c3*x*x*x        #Third order approximation
            if psitype == 'mtr4':
                psi_vals += (c3+c4*x)*x*x*x #Fourth order approximation

            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center] = finv(T_hat,ker,dx,F)
            if norm:T_vals[center] *= nrm(r,w_func,F,error_check=False)
            if verbose:
                print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Inversion"\
                % (i,n_used,nw,loop),end="\r")
    elif psitype == 'full':
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            d      = D_vals[center]
            b      = B_vals[center]
            phi0   = phi_rad_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            F      = F_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx,error_check=False)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                dphi_s_rad = psi_factor(r,r0,b,phi0,error_check=False)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
            else:
                crange    += 1
                r          = rho_vals[crange]
                phi_s_rad  = phi_s_rad1
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
            loop = 0
            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                loop      += 1
                if loop > 5:
                    break
            phi_s_rad1 = phi_s_rad

            # Compute psi and then compute averages and differences across psi.
            psi_vals = kD * psi(r,r0,d,b,phi_s_rad,phi0,error_check=True)
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center] = finv(T_hat,ker,dx,F)
            if norm:
                T_vals[center] *= nrm(r,psi_vals,w_func,F,error_check=False)
            if verbose:
                print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Inversion"\
                % (i,n_used,nw,loop),end="\r")
    else: raise TypeError("Illegal psitype: %s" % psitype)
    return T_vals

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

def fresnel_transform(T,ker,DX,f_scale):
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

def fresnel_inverse(T_hat,ker,dx,f_scale,error_check=True):
    """
        Function: fresnel_transform
        Purpose:  Compute the approximate inverse of a Fresnel transform.
            Variables:
                T_hat:   The input function (Diffraction data).
                KER:     The Fresnel kernel.
                dx:      The size of a bin (r[1]-r[0]).
                F_SCALE: The Frensel scale. If not set, a default of 1 is used.
        Output:
            T: The Fresnel inverted complex transmittance.
        History:
            Translated from IDL: RJM - 2018/04/15 12:23 P.M.
    """
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
    T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
    return T

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

def psi_factor(r,r0,b,phi0,error_check=True):
    """
        Function: psi_factor_fast
        Purpose:  Calculate the first iteration of Newton-Raphson for psi with
            respect to phi using previously calculated sines and cosines.
        Variables:
            r:    Ring radius variable.
            r0:   Ring intercept point.
            b:    Rring opening angle.
            phi0: Ring azimuth angle.
        History:
            Translated from IDL: RJM - Rough Draft - 2018/05/15 7:38 P.M.
    """
    if error_check:
        if (not check_real(r)):
            sys.exit("r must be real valued")
        if (not check_real(r0)):
            sys.exit("r0 must be real valued")
        if (not check_real(b)):
            sys.exit("b must be real valued")
        if (not check_real(phi0)):
            sys.exit("phi0 must be real valued")
    else: pass
    cb      = np.cos(b)
    sp0     = np.sin(phi0)
    cp0     = np.cos(phi0)
    factor  = ((cb*cb) * cp0 * sp0 / (1.0 - (cb*cb) * (sp0*sp0))) * (r - r0) / r0
    return factor

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