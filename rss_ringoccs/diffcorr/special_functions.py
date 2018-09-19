import numpy as np
from scipy.special import erf, lambertw
SQRT_PI_2 = 0.886226925452758013649084

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
        Smooth (and optionally differentiate) data with a
        Savitzky-Golay filter. The Savitzky-Golay filter removes
        high frequency noise from data. It has the advantage of
        preserving the original shape and features of the signal
        better than other types of filtering approaches, such as
        moving averages techniques.
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

    order_range = range(order+1)
    half_window = (window_size -1) // 2

    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = 1
    for i in range(deriv):
        m *= rate*(deriv-i)
    m *= np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

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

def sq_well_solve(x,a,b,F,invert=False):
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
