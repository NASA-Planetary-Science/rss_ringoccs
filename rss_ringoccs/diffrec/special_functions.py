import numpy as np
from scipy.special import erf, lambertw
from . import window_functions
from rss_ringoccs.tools import error_check
try:
    from rss_ringoccs._ufuncs import _special_functions
except (ImportError, ModuleNotFoundError):
    print(
        """
            Error: rss_ringoccs.diffrec.special_functions
            \tCould Not Import C Code. Stricly Using Python Code.
            \tThis is signicantly slower. There was most likely an error
            \tin your installation of rss_ringoccs. To use the C Code,
            \tdownload a C Compiler (GCC) and see the User's Guide for
            \tinstallation instructions.
        """
    )

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
        Purpose:
            To smooth data with a Savitzky-Golay filter.
            This removes high frequency noise while
            maintaining many of the original features of
            the input data.
        Arguments:
            :y (*np.ndarray*):
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
            :y_smooth (*np.ndarray*):
                The data smoothed by the Savitzky-Golay filter.
                This returns the nth derivative if the deriv
                keyword has been set.
        Notes:
            The Savitzky-Golay is a type of low-pass filter,
            particularly suited for smoothing noisy data.
            The main idea behind this approach is to make for
            each point a least-square fit with a polynomial of
            high order over a odd-sized window centered at the point.
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

def compute_norm_eq(w_func, error_check=True):
    """
        Purpose:
            Compute normalized equivalenth width of a given function.
        Arguments:
            :w_func (*np.ndarray*):
                Function with which to compute
                the normalized equivalent width of.
        Outputs:
            :normeq (*float*):
                The normalized equivalent width of w_func.
        Notes:
            The normalized equivalent width is effectively computed
            using Riemann sums to approximate integrals. Therefore
            large dx values (Spacing between points in w_func)
            will result in an inaccurate normeq. One should keep
            this in mind during calculations.
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

def fresnel_scale(Lambda, d, phi, b, deg=False):
    """
        Purpose:
            Compute the Fresnel Scale from lambda, D, Phi, and B.
        Arguments:
            :Lambda (*np.ndarray* or *float*):
                Wavelength of the incoming signal.
            :d (*np.ndarray* or *float*):
                RIP-Spacecraft Distance.
            :phi (*np.ndarray* or *float*):
                Ring azimuth angle.
            :b (*np.ndarray* or *float*):
                Ring opening angle.
        Keywords:
            :deg (*bool*):
                Set True if phi/b are in degrees.
                Default is radians.
        Output:
            :fres (*np.ndarray* or *float*):
                The Fresnel scale.
        Note:
            Lambda and d must be in the same units.
            The output (Fresnel scale) will have the same units
            as lambda and d. In addition, b and phi must also
            have the same units. If b and phi are in degrees,
            make sure to set deg=True. Default is radians.
    """
    try:
        Lambda = np.array(Lambda)
        phi = np.array(phi)
        d = np.array(d)
        b = np.array(b)

        if deg:
            b = np.deg2rad(b)
            phi = np.deg2rad(phi)
            cb = np.cos(b)
            sb = np.sin(b)
            sp = np.sin(phi)
        else:
            cb = np.cos(b)
            sb = np.sin(b)
            sp = np.sin(phi)

        return np.sqrt(0.5 * Lambda * d * (1-np.square(cb*sp)) / np.square(sb))
    except (TypeError, ValueError):
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\tdiffrec.special_functions.fresnel_scale\n
                \r\tInputs should be four numpy arrays or
                \r\tfloating point numbers.
            """
        )

def psi(kD, r, r0, phi, phi0, B, D):
    """
        Purpose:
            Compute psi (MTR Equation 4)
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*np.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*np.ndarray*):
                Root values of dpsi/dphi, radians.
            :phi0 (*np.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :psi (*np.ndarray*):
                Geometric Function from Fresnel Kernel.
    """
    # Compute Xi variable (MTR86 Equation 4b). Signs of xi are swapped.
    xi = (np.cos(B)/D) * (r * np.cos(phi) - r0 * np.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*np.cos(phi-phi0)) / (D*D)

    # Sign of xi swapped from MTR86.
    psi_vals = kD * (np.sqrt(1.0+eta-2.0*xi) + xi - 1.0)
    return psi_vals

def dpsi(kD, r, r0, phi, phi0, B, D):
    """
        Purpose:
            Compute dpsi/dphi
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*np.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*np.ndarray*):
                Root values of dpsi/dphi, radians.
            :phi0 (*np.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :dpsi (*array*):
                Partial derivative of psi with
                respect to phi.
    """
    # Compute Xi variable (MTR86 Equation 4b).
    xi = (np.cos(B)/D) * (r * np.cos(phi) - r0 * np.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*np.cos(phi-phi0)) / (D*D)

    psi0 = np.sqrt(1.0+eta-2.0*xi)

    # Compute derivatives.
    dxi = -(np.cos(B)/D) * (r*np.sin(phi))
    deta = 2.0*r*r0*np.sin(phi-phi0)/(D*D)

    # Compute the partial derivative.
    psi_d1 = (0.5/psi0)*(deta-2.0*dxi) + dxi
    psi_d1 *= kD

    return psi_d1

def dpsi_ellipse(kD, r, r0, phi, phi0, B, D, ecc, peri):
    """
        Purpose:
            Compute dpsi/dphi
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*np.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*np.ndarray*):
                Root values of dpsi/dphi, radians.
            :phi0 (*np.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :dpsi (*array*):
                Partial derivative of psi with
                respect to phi.
    """
    # Compute Xi variable (MTR86 Equation 4b).
    xi = (np.cos(B)/D) * (r * np.cos(phi) - r0 * np.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*np.cos(phi-phi0)) / (D*D)

    psi0 = np.sqrt(1.0+eta-2.0*xi)

    # Compute derivatives.
    dxi_phi = -(np.cos(B)/D) * (r*np.sin(phi))
    deta_phi = 2.0*r*r0*np.sin(phi-phi0)/(D*D)

    dxi_rho = (np.cos(B)/D)*np.cos(phi)
    deta_rho = 2.0*(r-r0*np.cos(phi-phi0)) / (D*D)

    # Compute the partial derivative.
    psi_d1 = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho
    psi_d1 *= r*ecc*np.sin(phi-peri)/(1+ecc*np.cos(phi-peri))
    psi_d1 += (deta_phi-2.0*dxi_phi)*(0.5/psi0) + dxi_phi

    psi_d1 *= kD

    return psi_d1

def d2psi(kD, r, r0, phi, phi0, B, D):
    """
        Purpose:
            Compute d^2psi/dphi^2
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*np.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*np.ndarray*):
                Root values of dpsi/dphi, radians.
            :phi0 (*np.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :dpsi (*np.ndarray*):
                Second partial derivative of psi
                with respect to phi.
    """
    # Compute Xi variable (MTR86 Equation 4b).
    xi = (np.cos(B)/D) * (r * np.cos(phi) - r0 * np.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*np.cos(phi-phi0)) / (D*D)

    psi0 = np.sqrt(1.0+eta-2.0*xi)

    # Compute derivatives.
    dxi = -(np.cos(B)/D) * (r*np.sin(phi))
    dxi2 = -(np.cos(B)/D) * (r*np.cos(phi))

    deta = 2.0*r*r0*np.sin(phi-phi0)/(D*D)
    deta2 = 2.0*r*r0*np.cos(phi-phi0)/(D*D)

    # Compute the second partial derivative.
    psi_d2 = (-0.25/(psi0*psi0*psi0))*(deta-2.0*dxi)*(deta-2.0*dxi)
    psi_d2 += (0.5/psi0)*(deta2-2.0*dxi2)+dxi2
    psi_d2 *= kD

    return psi_d2

def resolution_inverse(x):
    """
        Purpose:
            Compute the inverse of y = x/(exp(-x)+x-1)
        Arguments:
            :x (*np.ndarray* or *float*):
                Independent variable
        Outputs:
            :f (*np.ndarray* or *float*):
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
    fname = "diffrec.special_functions.resolution_inverse"
    error_check.check_is_real(x, "x", fname)
    error_check.check_positive(np.min(np.real(x))-1.0, "x-1", fname)

    P1 = x/(1.0-x)
    P2 = P1*np.exp(P1)
    f = lambertw(P2)-P1

    if np.all(np.isreal(x)):
        f = np.real(f)

    return f

def fresnel_cos(x):
    """
        Purpose:
            Compute the Fresnel cosine function.
        Arguments:
            :x (*np.ndarray* or *float*):
                A real or complex number, or numpy array.
        Outputs:
            :f_cos (*np.ndarray* or *float*):
                The fresnel cosine integral of x.
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
        Examples:
            Compute and plot the Fresnel Cosine integral.
                In [1]: import rss_ringoccs.diffcorr.special_functions as sf
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = sf.fresnel_cos(x)
                In [6]: plt.show(plt.plot(x,y))
    """
    try:
        return _special_functions.fresnel_cos(x)
    except (TypeError, ValueError, NameError):
        fname = "diffrec.special_functions.fresnel_cos"
        error_check.check_is_real(x, "x", fname)

        x *= window_functions.RCP_SQRT_2
        f_cos = ((0.25-0.25j)*erf((1.0+1.0j)*x)+(0.25+0.25j)*erf((1.0-1.0j)*x))

        if (np.isreal(x).all()):
            f_cos = np.real(f_cos)

        return f_cos*window_functions.SQRT_PI_2

def fresnel_sin(x):
    """
        Purpose:
            Compute the Fresnel sine function.
        Variables:
            :x (*np.ndarray* or *float*):
                The independent variable.
        Outputs:
            :f_sin (*np.ndarray* or *float*):
                The fresnel sine integral of x.
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
        Examples:
            Compute and plot the Fresnel Sine integral.
                In [1]: import rss_ringoccs.diffcorr.special_functions as sf
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = sf.fresnel_sin(x)
                In [6]: plt.show(plt.plot(x,y))
    """
    try:
        return _special_functions.fresnel_sin(x)
    except (TypeError, ValueError, NameError):
        fname = "diffrec.special_functions.fresnel_sin"
        error_check.check_is_real(x, "x", fname)

        x *= window_functions.RCP_SQRT_2
        f_sin = ((0.25+0.25j)*erf((1.0+1.0j)*x)+(0.25-0.25j)*erf((1.0-1.0j)*x))

        if (np.isreal(x).all()):
            f_sin = np.real(f_sin)

        return f_sin*window_functions.SQRT_PI_2

def square_well_diffraction(x, a, b, F):
    try:
        return _special_functions.square_well_diffraction(x, a, b, F)
    except(TypeError, ValueError, NameError):
        try:
            arg_1 = np.sqrt(np.pi/2.0)*((a-x)/F)
            arg_2 = np.sqrt(np.pi/2.0)*((b-x)/F)
            
            return 1.0 - np.sqrt(2.0/np.pi)*(0.5 - 0.5j) * (
                fresnel_cos(arg_2)-fresnel_cos(arg_1)+
                1j*(fresnel_sin(arg_2)-fresnel_sin(arg_1))
            )
        except(TypeError, ValueError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\tdiffrec.special_functions.square_well_diffraction
                    \r
                    \r\tInvalid input. Input should be:
                    \r\t\tx:\t Numpy array of floating point numbers.
                    \r\t\ta:\t Floating point number
                    \r\t\tb:\t Floating point number
                    \r\t\tF:\t Floating point number
                """
            )

def inverse_square_well_diffraction(x, a, b, F):
    try:
        return _special_functions.inverse_square_well_diffraction(x, a, b, F)
    except(TypeError, ValueError, NameError):
        try:
            arg_1 = np.sqrt(np.pi/2.0)*((a-x)/F)
            arg_2 = np.sqrt(np.pi/2.0)*((b-x)/F)
            
            return np.sqrt(2.0/np.pi)*(0.5 - 0.5j) * (
                fresnel_cos(arg_2)-fresnel_cos(arg_1)+
                1j*(fresnel_sin(arg_2)-fresnel_sin(arg_1))
            )
        except(TypeError, ValueError):
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\tdiffrec.special_functions.inverse_square_well_diffraction
                    \r
                    \r\tInvalid input. Input should be:
                    \r\t\tx:\t Numpy array of floating point numbers.
                    \r\t\ta:\t Floating point number
                    \r\t\tb:\t Floating point number
                    \r\t\tF:\t Floating point number
                """
            )

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
                The perpendicular distance from the slit plane to
                the observer.
            :a (*float*):
                The slit parameter. This is a unitless paramter
                defined as the ratio between the slit width and
                the wavelength of the incoming signal.
        Outputs:
            :f:
                Single slit diffraction pattern.
    """
    fname = "diffrec.special_functions.single_slit_diffraction"
    error_check.check_is_real(x, "x", fname)
    z = error_check.check_type_and_convert(z, float, "z", fname)
    a = error_check.check_type_and_convert(a, float, "a", fname)

    return np.square(np.sinc(a*x/z))

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
    fname = "diffrec.special_functions.double_slit_diffraction"
    error_check.check_is_real(x, "x", fname)
    z = error_check.check_type_and_convert(z, float, "z", fname)
    a = error_check.check_type_and_convert(a, float, "a", fname)
    f1 = np.square(np.sinc(a*x/z))
    f2 = np.square(np.sin(window_functions.TWO_PI*d*x/z))
    f3 = 4.0*np.square(np.sin(np.pi*d*x/z))

    return f1*f2/f3
