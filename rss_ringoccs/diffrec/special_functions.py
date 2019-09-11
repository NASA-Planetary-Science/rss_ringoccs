import numpy
from . import window_functions
from rss_ringoccs.tools import error_check
try:
    from . import _special_functions, _diffraction_functions
except:
    raise ImportError(
        """
        \r\tError: rss_ringoccs
        \r\t\tdiffrec.special_functions\n
        \r\tCould Not Import C Code. There was most likely an error
        \r\tin your installation of rss_ringoccs. Install GCC (C Compiler)
        \r\tand see the User's Guide for installation instructions.
        """
    )

def besselJ0(x):
    """
        Purpose:
            Compute the Bessel Function J_0(x).
        Variables:
            :x (*numpy.ndarray* or *float*):
                The independent variable.
        Outputs:
            :J_0 (*numpy.ndarray* or *float*):
                The fresnel sine integral of x.
        Notes:

        Examples:
            Compute the Fresnel Sine integral from -10 to 10, with
            points spaced 0.01 apart.

            >>> import rss_ringoccs.diffcorr.special_functions as sf
            >>> import numpy as np
            >>> x = numpy.arange(-10, 10, 0.01)
            >>> y = sf.besselJ0(x)
    """
    try:
        return _special_functions.besselJ0(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.besselJ0\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = besselJ0(x)
            """
        )

def besselI0(x):
    """
        Purpose:
            Compute the Bessel Function J_0(x).
        Variables:
            :x (*numpy.ndarray* or *float*):
                The independent variable.
        Outputs:
            :J_0 (*numpy.ndarray* or *float*):
                The fresnel sine integral of x.
        Notes:

        Examples:
            Compute the Fresnel Sine integral from -10 to 10, with
            points spaced 0.01 apart.

            >>> import rss_ringoccs.diffcorr.special_functions as sf
            >>> import numpy as np
            >>> x = numpy.arange(-10, 10, 0.01)
            >>> y = sf.besselI0(x)
    """
    try:
        return _special_functions.besselI0(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.besselI0\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = besselI0(x)
            """
        )

def frequency_to_wavelength(freq_hz):
    try:
        return _special_functions.frequency_to_wavelength(freq_hz)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.wavelength_to_wavenumber\n
            \r\tInput should be a numpy array of non-zero real numbers
            \r\t(ints or floats), or a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(3, 10)
            \r\t\t>>> y = frequency_to_wavelength(x)
            """
        )

#TODO
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
    Purpose:
        To smooth data with a Savitzky-Golay filter. This removes
        high frequency noise while maintaining many of the
        original features of the input data.
    Arguments:
        :y (*numpy.ndarray*):
            The input "Noisy" data.
        :window_size (*int*):
            The length of the window. Must be an odd number.
        :order (*int*):
            The order of the polynomial used for filtering.
            Must be less then window_size - 1.
    Keywords:
        :deriv (*int*):
            The order of the derivative what will be computed.
    Output:
        :y_smooth (*numpy.ndarray*):
            The data smoothed by Savitzky-Golay filter.
    """
    try:
        y = numpy.array(y)
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs
            \r\t\tdiffrec.special_functions.savitzky_golay\n
            \r\tUsage:
            \r\t\tysmooth = savitzky_golay(y, Window_Size, Poly_Order)
            \r\t\ty:            Numpy Array (Floating Point)
            \r\t\tWindow_Size   Positive Odd Integer
            \r\t\tPoly_Order    Positive Integer (Poly_Order > Window_Size)
            """
        )

    if (window_size % 2 != 1) or (window_size < 1):
        raise ValueError(
            """
            \r\tError Encountered: rss_ringoccs
            \r\t\tdiffrec.special_functions.savitzky_golay\n
            \r\twindow_size must be an odd integer.
            """
        )
    elif (window_size < order + 2):
        raise ValueError(
            """
            \r\tError Encountered: rss_ringoccs
            \r\t\tdiffrec.special_functions.savitzky_golay\n
            \r\twindow_size must be less than Poly_Order.
            """
        )
    else:
        pass

    half_window = (window_size - 1) // 2

    # precompute coefficients
    b = numpy.zeros((window_size, order+1))
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

    b = numpy.mat(b)

    m = 1
    for i in range(deriv):
        m *= rate*(deriv-i)

    m *= numpy.linalg.pinv(b).A[deriv]

    # Pad the endpoints with values from the signal.
    firstvals = y[0] - numpy.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))

    return numpy.convolve(m[::-1], y, mode='valid')

def compute_norm_eq(w_func):
    """
    Purpose:
        Compute normalized equivalenth width of a given function.
    Arguments:
        :w_func (*numpy.ndarray*):
            Function to compute the normalized equivalent width.
    Outputs:
        :normeq (*float*):
            The normalized equivalent width of w_func.
    Notes:
        The normalized equivalent width is computed using Riemann
        sums to approximate integrals. Therefore large dx values
        (Spacing between points) will result in an inaccurate
        normeq. One should keep this in mind during calculations.
    Examples:
        Compute the Kaiser-Bessel 2.5 window of width 30km and
        spacing 0.1 and compute the normalized equivalent width:
            >>> from rss_ringoccs import diffrec as dc
            >>> w = dc.window_functions.kb25(30, 0.1)
            >>> normeq = dc.special_functions.compute_norm_eq(w)
            >>> print(normeq)
            1.6573619266424229
        In contrast, the actual value is 1.6519208.
        Compute the normalized equivalent width for the squared
        cosine window of width 10 and spacing 0.25.
            >>> from rss_ringoccs import diffrec as dc
            >>> w = dc.window_functions.coss(10, 0.25)
            >>> normeq = dc.special_functions.compute_norm_eq(w)
            >>> print(normeq)
            1.5375000000000003
        The normalized equivalent width of the squared cosine
        function can be computed exactly using standard methods
        from a calculus course. It's value is exactly 1.5
        If we use a smaller dx when computing w, we get a better
        approximation. Use width 10 and spacing 0.001.
            >>> from rss_ringoccs import diffrec as dc
            >>> w = dc.window_functions.coss(10, 0.001)
            >>> normeq = dc.special_functions.compute_norm_eq(w)
            >>> print(normeq)
            1.50015
    """
    try:
        return _special_functions.compute_norm_eq(w_func)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.wavelength_to_wavenumber\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = frequency_to_wavelength(x)
            """
        )

def fresnel_scale(Lambda, d, phi, b):
    """
    Purpose:
        Compute the Fresnel Scale.
    Arguments:
        :Lambda (*numpy.ndarray* or *float*):
            Wavelength of the incoming signal.
        :d (*numpy.ndarray* or *float*):
            RIP-Spacecraft Distance.
        :phi (*numpy.ndarray* or *float*):
            Ring azimuth angle.
        :b (*numpy.ndarray* or *float*):
            Ring opening angle.
    Output:
        :fres (*numpy.ndarray* or *float*):
            The Fresnel scale.
    Note:
        :math:`\\lambda` and :math:`D` must be in the same units.
        The output (Fresnel scale) will have the same units
        as :math:`\\lambda` and d. In addition, :math:`B`
        and :math:`\\phi` must also be in radians.
    """
    try:
        return _special_functions.fresnel_scale(Lambda, d, phi, b)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_scale\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> Lambda = 3.6e-6
            \r\t\t>>> d = 300000.0
            \r\t\t>>> phi = numpy.arange(0, 3.14, 0.01)
            \r\t\t>>> b = 0.3
            \r\t\t>>> y = fresnel_scale(Lambda, d, phi, b)\n
            """
        )

def fresnel_psi(kD, r, r0, phi, phi0, B, D):
    """
    Purpose:
        Compute :math:`\\psi` (MTR Equation 4)
    Arguments:
        :kD (*float*):
            Wavenumber, unitless.
        :r (*float*):
            Radius of reconstructed point, in kilometers.
        :r0 (*numpy.ndarray*):
            Radius of region within window, in kilometers.
        :phi (*numpy.ndarray*):
            Ring azimuth angle corresponding to r, radians.
        :phi0 (*numpy.ndarray*):
            Ring azimuth angle corresponding to r0, radians.
        :B (*float*):
            Ring opening angle, in radians.
        :D (*float*):
            Spacecraft-RIP distance, in kilometers.
    Outputs:
        :psi (*numpy.ndarray*):
            Geometric Function from Fresnel Kernel.
    """
    try:
        return _special_functions.fresnel_psi(kD, r, r0, phi, phi0, B, D)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_psi\n
            \r\tInput should be seven numpy arrays of real numbers.\n
            \r\tUsage:
            \r\t\t>>> y = fresnel_psi(kD, r, r0, phi, phi0, B, D)\n
            """
        )

def fresnel_dpsi_dphi(kD, r, r0, phi, phi0, B, D):
    """
        Purpose:
            Compute :math:`\\mathrm{d}\\psi/\\mathrm{d}\\phi`
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*numpy.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*numpy.ndarray*):
                Ring azimuth angle corresponding to r, radians.
            :phi0 (*numpy.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :dpsi (*array*):
                Partial derivative of :math:`\\psi` with
                respect to :math:`\\phi`.
    """
    try:
        return _special_functions.fresnel_dpsi_dphi(kD, r, r0, phi, phi0, B, D)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_dpsi_dphi\n
            \r\tInput should be seven numpy arrays of real numbers.\n
            \r\tUsage:
            \r\t\t>>> y = fresnel_dpsi_dphi(kD, r, r0, phi, phi0, B, D)\n
            """
        )

# TODO
def dpsi_ellipse(kD, r, r0, phi, phi0, B, D, ecc, peri):
    """
        Purpose:
            Compute :math:`\\mathrm{d}\\psi/\\mathrm{d}\\phi`
        Arguments:
            :kD (*float*):
                Wavenumber, unitless.
            :r (*float*):
                Radius of reconstructed point, in kilometers.
            :r0 (*numpy.ndarray*):
                Radius of region within window, in kilometers.
            :phi (*numpy.ndarray*):
                Root values of :math:`\\mathrm{d}\\psi/\\mathrm{d}\\phi`, radians.
            :phi0 (*numpy.ndarray*):
                Ring azimuth angle corresponding to r0, radians.
            :B (*float*):
                Ring opening angle, in radians.
            :D (*float*):
                Spacecraft-RIP distance, in kilometers.
        Outputs:
            :dpsi (*array*):
                Partial derivative of :math:`\\psi` with
                respect to :math:`\\phi`.
    """
    # Compute Xi variable (MTR86 Equation 4b).
    xi = (numpy.cos(B)/D) * (r * numpy.cos(phi) - r0 * numpy.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*numpy.cos(phi-phi0)) / (D*D)

    psi0 = numpy.sqrt(1.0+eta-2.0*xi)

    # Compute derivatives.
    dxi_phi = -(numpy.cos(B)/D) * (r*numpy.sin(phi))
    deta_phi = 2.0*r*r0*numpy.sin(phi-phi0)/(D*D)

    dxi_rho = (numpy.cos(B)/D)*numpy.cos(phi)
    deta_rho = 2.0*(r-r0*numpy.cos(phi-phi0)) / (D*D)

    # Compute the partial derivative.
    psi_d1 = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho
    psi_d1 *= r*ecc*numpy.sin(phi-peri)/(1+ecc*numpy.cos(phi-peri))
    psi_d1 += (deta_phi-2.0*dxi_phi)*(0.5/psi0) + dxi_phi

    psi_d1 *= kD

    return psi_d1

# TODO
def fresnel_d2psi_dphi2(kD, r, r0, phi, phi0, B, D):
    """
    Purpose:
        Compute :math:`\\mathrm{d}^2\\psi/\\mathrm{d}\\phi^2`
    Arguments:
        :kD (*float*):
            Wavenumber, unitless.
        :r (*float*):
            Radius of reconstructed point, in kilometers.
        :r0 (*numpy.ndarray*):
            Radius of region within window, in kilometers.
        :phi (*numpy.ndarray*):
            Root values of :math:`\\mathrm{d}\\psi/\\mathrm{d}\\phi`,
            radians.
        :phi0 (*numpy.ndarray*):
            Ring azimuth angle corresponding to r0, radians.
        :B (*float*):
            Ring opening angle, in radians.
        :D (*float*):
            Spacecraft-RIP distance, in kilometers.
    Outputs:
        :dpsi (*numpy.ndarray*):
            Second partial derivative of :math:`\\psi`
            with respect to :math:`\\phi`.
    """
    # Compute Xi variable (MTR86 Equation 4b).
    xi = (numpy.cos(B)/D) * (r * numpy.cos(phi) - r0 * numpy.cos(phi0))

    # Compute Eta variable (MTR86 Equation 4c).
    eta = (r0*r0 + r*r - 2.0*r*r0*numpy.cos(phi-phi0)) / (D*D)

    psi0 = numpy.sqrt(1.0+eta-2.0*xi)

    # Compute derivatives.
    dxi = -(numpy.cos(B)/D) * (r*numpy.sin(phi))
    dxi2 = -(numpy.cos(B)/D) * (r*numpy.cos(phi))

    deta = 2.0*r*r0*numpy.sin(phi-phi0)/(D*D)
    deta2 = 2.0*r*r0*numpy.cos(phi-phi0)/(D*D)

    # Compute the second partial derivative.
    psi_d2 = (-0.25/(psi0*psi0*psi0))*(deta-2.0*dxi)*(deta-2.0*dxi)
    psi_d2 += (0.5/psi0)*(deta2-2.0*dxi2)+dxi2
    psi_d2 *= kD

    return psi_d2

def resolution_inverse(x):
    """
    Purpose:
        Compute the inverse of :math:`y = x/(\\exp(-x)+x-1)`
    Arguments:
        :x (*numpy.ndarray* or *float*):
            Independent variable
    Outputs:
        :f (*numpy.ndarray* or *float*):
            The inverse of :math:`x/(\\exp(-x)+x-1)`
    Dependencies:
        #. numpy
        #. scipy.special
    Method:
        The inverse of :math:`x/(\\exp(-x)+x-1)` is computed using the
        LambertW function. This function is the inverse of
        :math:`y = x\\exp(x)`. This is computed using the scipy.special
        subpackage using their lambertw function.
    Warnings:
        #. The real part of the argument must be greater than 1.
        #. The scipy.special lambertw function is slightly
           inaccurate when it's argument is near :math:`-1/e`. This
           argument is :math:`z = \\exp(x/(1-x)) * x/(1-x)`
    Examples:
        Plot the function on the interval (1,2)

        >>> import rss_ringoccs.diffcorr.special_functions as sf
        >>> import numpy as np
        >>> x = numpy.array(range(0,1001))*0.001+1.01
        >>> y = sf.resolution_inverse(x)
        >>> import matplotlib.pyplot as plt
        >>> plt.show(plt.plot(x,y))
    """
    try:
        return _special_functions.resolution_inverse(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.resolution_inverse\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor an int or a float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = resolution_inverse(x)
            """
        )

def square_well_diffraction(x, a, b, F):
    try:
        return _special_functions.square_well_diffraction(x, a, b, F)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.square_well_diffraction\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand three floats/ints.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-10, 10, 0.01)
            \r\t\t>>> a = -5.0
            \r\t\t>>> b = 5.0
            \r\t\t>>> F = 0.05
            \r\t\t>>> y = square_well_diffraction(x, a, b, F)
            """
        )

def inverse_square_well_diffraction(x, a, b, F):
    try:
        return _special_functions.inverse_square_well_diffraction(x, a, b, F)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.inverse_square_well_diffraction\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand three floats/ints.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-10, 10, 0.01)
            \r\t\t>>> a = -5.0
            \r\t\t>>> b = 5.0
            \r\t\t>>> F = 0.05
            \r\t\t>>> y = inverse_square_well_diffraction(x, a, b, F)
            """
        )

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
    """
    try:
        return _special_functions.double_slit_diffraction(x, z, a, d)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.double_slit_diffraction\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand three floats/ints.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-10, 10, 0.01)
            \r\t\t>>> z = 5.0
            \r\t\t>>> a = 10.0
            \r\t\t>>> d = 1.0
            \r\t\t>>> y = double_slit_diffraction(x, z, a, d)
            """
        )

def fresnel_cos(x):
    """
    Purpose:
        Compute the Fresnel cosine function.
    Arguments:
        :x (*numpy.ndarray* or *float*):
            A real or complex number, or numpy array.
    Outputs:
        :f_cos (*numpy.ndarray* or *float*):
            The fresnel cosine integral of x.
    Notes:
        #.  The Fresnel Cosine integral is the solution to the equation
            :math:`\\mathrm{d}y/\\mathrm{d}x = \\cos(\\frac\\pi 2 x^2)`,
            :math:`y(0) = 0`. In other words,
            :math:`y = \\int_{t=0}^{x}\\cos(\\frac\\pi 2 t^2)\\mathrm{d}t`
        #.  The Fresnel Cosine integral is used for the solution
            of diffraction through a square well. Because of this
            it is useful for forward modeling problems in
            radiative transfer and diffraction.
    Examples:
        Compute and the Fresnel Cosine integral from -10 to 10, with
            points spaced 0.01 apart.

        >>> import rss_ringoccs.diffcorr.special_functions as sf
        >>> import numpy as np
        >>> x = numpy.arange(-10, 10, 0.01)
        >>> y = sf.fresnel_cos(x)
    """
    try:
        return _special_functions.fresnel_cos(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_cos\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = fresnel_cos(x)
            """
        )

def fresnel_sin(x):
    """
        Purpose:
            Compute the Fresnel sine function.
        Variables:
            :x (*numpy.ndarray* or *float*):
                The independent variable.
        Outputs:
            :f_sin (*numpy.ndarray* or *float*):
                The fresnel sine integral of x.
        Notes:
            #.  The Fresnel sine integral is the solution to the equation
                :math:`\\mathrm{d}y/\\mathrm{d}x = \\sin(\\frac\\pi 2 x^2)`,
                :math:`y(0) = 0`. In other words,
                :math:`y = \\int_{t=0}^{x}\\sin(\\frac\\pi 2 t^2) dt`
            #.  The Fresnel sine integral is used for the solution
                of diffraction through a square well. Because of this
                is is useful for forward modeling problems in
                radiative transfer and diffraction.
        Examples:
            Compute the Fresnel Sine integral from -10 to 10, with
            points spaced 0.01 apart.

            >>> import rss_ringoccs.diffcorr.special_functions as sf
            >>> import numpy as np
            >>> x = numpy.arange(-10, 10, 0.01)
            >>> y = sf.fresnel_sin(x)
    """
    try:
        return _special_functions.fresnel_sin(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_sin\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = fresnel_sin(x)
            """
        )

def fresnel_transform(T_in, rho_km_vals, F_km_vals, w_km_vals, start, n_used,
                      wtype, norm, fwd, psitype, phi_rad_vals=None,
                      kD_vals=None, B_rad_vals=None, D_km_vals=None,
                      periapse=None, eccentricity=None):

    fname = "diffrec.special_functions.fresnel_transform"
    # Remove spaces/quotes from the wtype variable and set to lower case.
    wtype = wtype.replace(" ", "").replace("'", "").replace('"', "")
    wtype = wtype.lower()

    # Check that wtype is in the allowed list of window types.
    if not (wtype in window_functions.func_dict):
        erm = ""
        for key in window_functions.func_dict:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise ValueError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \tIllegal string used for wtype.
                \r\t\tYour string: '%s'\n
                \r\tAllowed Strings:\n%s
            """ % (fname, wtype, erm)
        )
    else:
        pass

    # Check that range and psitype are legal inputs.
    psitype = error_check.check_psitype(psitype, fname)

    if (psitype == "ellipse"):

        # Compute the sample spacing.
        dx_km = rho_km_vals[1]-rho_km_vals[0]

        # Extract the window function from the window function dictionary.
        fw = window_functions.func_dict[wtype]["func"]

        # Create empty array for reconstruction / forward transform.
        T_out = T_in * 0.0

        # Compute first window width and window function.
        w_init = w_km_vals[start]

        # Number of points in the first window (Odd integer).
        nw = int(2 * numpy.floor(w_init / (2.0 * dx_km)) + 1)

        # Indices for the data corresponding to current window.
        crange = numpy.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))

        # Various geometry variables.
        r0 = rho_km_vals[crange]
        r = rho_km_vals[start]
        x = r-r0

        # Compute the first window function.
        w_func = fw(x, w_init)

        # Perform the Fresnel Transform, point by point, via Riemann sums.
        for i in numpy.arange(n_used):

            # Current point being computed.
            center = start+i

            # Current window width, Fresnel scale, and ring radius.
            w = w_km_vals[center]
            F = F_km_vals[center]
            r = rho_km_vals[center]

            # If window widths changes too much, recompute the window function.
            if (numpy.abs(w_init - w) >= 2.0 * dx_km):

                # Compute first window width and window function.
                w_init = w_km_vals[center]

                # Number of points in window (Odd integer).
                nw = int(2 * numpy.floor(w_init / (2.0 * dx_km)) + 1)

                # Indices for the data corresponding to this window.
                crange = numpy.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))

                # Ajdust ring radius by dx_km.
                r0 = rho_km_vals[crange]
                x = r-r0

                # Recompute the window function.
                w_func = fw(x, w_init)

            else:

                # If width hasn't changed much, increment index variable.
                crange += 1
                r0 = rho_km_vals[crange]

            # Various geometry variables for the current point.
            d = D_km_vals[crange]
            b = B_rad_vals[crange]
            kD = kD_vals[crange]
            phi = phi_rad_vals[crange]
            phi0 = phi_rad_vals[crange]

            # Compute Newton-Raphson perturbation
            psi_d1 = dpsi_ellipse(kD, r, r0, phi, phi0, b,
                                  d, eccentricity, periapse)
            loop = 0
            while (numpy.max(numpy.abs(psi_d1)) > 1.0e-4):
                psi_d1 = dpsi_ellipse(kD, r, r0, phi, phi0,
                                      b, d, eccentricity, periapse)
                psi_d2 = fresnel_d2psi_dphi2(kD, r, r0, phi, phi0, b, d)

                # Newton-Raphson
                phi += -(psi_d1 / psi_d2)

                # Add one to loop variable for each iteration
                loop += 1
                if (loop > 4):
                    break

            # Compute Psi (Fresnel Kernel, MTR86 Equation 4).
            psi_vals = fresnel_psi(kD, r, r0, phi, phi0, b, d)

            # Compute kernel function for Fresnel inverse or forward model.
            if fwd:
                ker = w_func*numpy.exp(1j*psi_vals)
            else:
                ker = w_func*numpy.exp(-1j*psi_vals)

            # Compute approximate Fresnel transform for current point.
            T = T_in[crange]
            T_out[center] = numpy.sum(ker*T) * dx_km * (0.5+0.5j)/F

            # If normalization has been set, normalize the reconstruction
            if norm:
                T_out[center] *= window_functions.window_norm(ker, dx_km, F)
        return T_out
    elif (psitype == "fresnel"):
        ord = 1
    elif (psitype == "fresnel4"):
        ord = 3
    elif (psitype == "fresnel6"):
        ord = 5
    elif (psitype == "fresnel8"):
        ord = 7
    elif (psitype == "full"):
        ord = 0
    else:
        raise ValueError(
            """
            \r\tError Encountered: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_transform\n
            \r\tIllegal psitype. Allowed options:
            \r\t\t'fresnel', 'fresnel4', 'fresnel6', 'fresnel8', 'full'
            """
        )

    # Try using the C version of the Fresnel transform.
    try:

        # Compute the Fresnel transform.
        return _diffraction_functions.fresnel_transform(
            T_in, rho_km_vals, F_km_vals, phi_rad_vals, kD_vals, B_rad_vals,
            D_km_vals, w_km_vals, start, n_used,
            window_functions.func_dict[wtype]["wnum"], int(norm), int(fwd), ord
        )
    except KeyboardInterrupt:
        raise
    except Exception as mes:
        print(mes)
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs
            \r\t\tdiffrec.special_functions.fresnel_transform\n
            \r\tArguments:
            \r\t\tT_in:         Complex numpy array.
            \r\t\trho_km_vals:  Positive real valued numpy array.
            \r\t\tF_km_vals:    Positive real valued numpy array.
            \r\t\tphi_rad_vals  Real valued numpy array.
            \r\t\tkD_vals:      Positive real valued numpy array.
            \r\t\tB_rad_vals:   Real valued numpy array.
            \r\t\tD_km_vals:    Positive real valued numpy array.
            \r\t\tw_km_vals:    Positive real valued numpy array.
            \r\t\tstart:        Integer.
            \r\t\tn_used:       Integer.
            \r\t\twtype:        String, name of the selected window.
            \r\t\tnorm:         Boolean.
            \r\t\tfwd:          Boolean.
            \r\t\tpsitype:      String.
            """
        )

def lambertw(x):
    try:
        return _special_functions.lambertw(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.lambertw\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor an int or a float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = lambertw(x)
            """
        )

def left_straightedge(x, edge, F):
    try:
        return _special_functions.left_straightedge(x, edge, F)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.left_straightedge\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand two positive real numbers.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(0, 1000, 0.1)
            \r\t\t>>> edge = 500.0
            \r\t\t>>> F = 2.0
            \r\t\t>>> y = left_straightedge(x, edge, F)
            """
        )

def max(x):
    """
    Purpose:
        Compute the maximum of a one dimensional numpy array.
        This function was written to test use of the C-Python API.
    Arguments:
        :x (*numpy.ndarray*):
            A one dimensional numpy array of real numbers.
    Outputs:
        :max (*float* or *int*):
            The maximum value of x.
    """
    try:
        return _special_functions.max(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.max\n
            \r\tInput should be a numpy array of real numbers.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.random.rand(100)
            \r\t\t>>> y = max(x)\n
            """
        )

def min(x):
    """
    Purpose:
        Compute the maximum of a one dimensional numpy array.
        This function was written to test use of the C-Python API.
    Arguments:
        :x (*numpy.ndarray*):
            A one dimensional numpy array of real numbers.
    Outputs:
        :max (*float* or *int*):
            The maximum value of x.
    """
    try:
        return _special_functions.min(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.min\n
            \r\tInput should be a numpy array of real numbers.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.random.rand(100)
            \r\t\t>>> y = min(x)\n
            """
        )

def right_straightedge(x, edge, F):
    try:
        return _special_functions.right_straightedge(x, edge, F)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.right_straightedge\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand two positive real numbers.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(0, 1000, 0.1)
            \r\t\t>>> edge = 500.0
            \r\t\t>>> F = 2.0
            \r\t\t>>> y = right_straightedge(x, edge, F)
            """
        ) 

def sinc(x):
    try:
        return _special_functions.sinc(x)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.sinc\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tor an int or a float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(-5, 5)
            \r\t\t>>> y = sinc(x)
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
    try:
        return _special_functions.single_slit_diffraction(x, z, a)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.single_slit_diffraction\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand three floats/ints.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-10, 10, 0.01)
            \r\t\t>>> z = 5.0
            \r\t\t>>> a = 10.0
            \r\t\t>>> y = single_slit_diffraction(x, z, a)
            """
        )

def square_well_phase(x, a, b, F):
    try:
        return _special_functions.square_well_phase(x, a, b, F)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.square_well_phase\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand three floats/ints.\n
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-10, 10, 0.01)
            \r\t\t>>> a = -5.0
            \r\t\t>>> b = 5.0
            \r\t\t>>> F = 0.05
            \r\t\t>>> y = square_well_phase(x, a, b, F)
            """
        )

def wavelength_to_wavenumber(lambda_km):
    try:
        return _special_functions.wavelength_to_wavenumber(lambda_km)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.special_functions.wavelength_to_wavenumber\n
            \r\tInput should be a numpy array of non-zero real numbers
            \r\t(ints or floats), or a non-zero int or non-zero float.\n
            \r\tUsage:
            \r\t\t>>> x = 1.0   # Or a numpy array, i.e. numpy.arange(3, 10)
            \r\t\t>>> y = wavelength_to_wavenumber(x)
            """
        )
