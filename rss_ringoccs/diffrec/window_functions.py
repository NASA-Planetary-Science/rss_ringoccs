"""
    Purpose:
        Provide a suite of window functions and functions related
        to the normalized equivalent width of a given array.
    Dependencies:
        #. numpy
        #. spicy
"""

import numpy as np
from scipy.special import lambertw, iv
from rss_ringoccs.tools import error_check
try:
    from . import _window_functions, special_functions
except:
    raise ImportError(
        """
        \r\tError: rss_ringoccs
        \r\t\tdiffrec.window_functions\n
        \r\tCould Not Import C Code. There was most likely an error
        \r\tin your installation of rss_ringoccs. Install GCC (C Compiler)
        \r\tand see the User's Guide for installation instructions.
        """
    )

# Declare constants for multiples of pi.
TWO_PI = 6.283185307179586476925287
SQRT_PI_2 = 1.253314137315500251207883
RCP_SQRT_2 = 0.7071067811865476

# Declare constants for the reciprocal of e and the square root of 2.
RCPR_E = 0.3678794411714423215955238
SQRT_2 = 1.414213562373095048801689

def rect(x, W):
    """
        Purpose:
            Compute the rectangular window function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.rect(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.rect\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = rect(x, W)
            """
        )

def coss(x, W):
    """
        Purpose:
            Compute the Cosine Squared Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.coss(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.coss\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = coss(x, W)
            """
        )

def kb20(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 2.0 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kb20(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kb20\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kb20(x, W)
            """
        )

def kb25(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 2.5 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kb25(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kb25\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kb25(x, W)
            """
        )

def kb35(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 3.5 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kb35(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kb35\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kb35(x, W)
            """
        )

def kbmd20(x, W):
    """
        Purpose:
            Compute the Modified Kaiser-Bessel 2.0 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kbmd20(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kbmd20\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kbmd20(x, W)
            """
        )

def kbmd25(x, W):
    """
        Purpose:
            Compute the Modified Kaiser-Bessel 2.5 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kbmd25(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kbmd25\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kbmd25(x, W)
            """
        )

def kbmd35(x, W):
    """
        Purpose:
            Compute the Modified Kaiser-Bessel 3.5 Window Function.
        Arguments:
            :x (*np.ndarray*):
                Real valued array used for the independent
                variable, or x-axis.
            :w_in (*float*):
                Width of the window. Positive float.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width w_in evaluated along x.
    """
    try:
        return _window_functions.kbmd35(x, W)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.kbmd35\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand a positive floating point number.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> W = 5.0
            \r\t\t>>> y = kbmd35(x, W)
            """
        )

def kbal(x, W, alpha):
    """
        Purpose:
            Create a Kaiser-Bessel window with a user
            specified alpha parameter.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
            :al (*float*):
                The alpha parameter :math:`\\alpha_0`.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel alpha window of width
                w_in and spacing dx between points.
        Notes:
            #.  The Kaiser-Bessel window is computed using the
                modified Bessel Function of the First Kind. It's
                value is :math:`y = I_0(\\alpha\\sqrt{1-4x^2/w^2})/I_0(\\alpha)`,
                where w is the window width.
            #.  We automatically multiply the alpha parameter by pi,
                so the kbal window function has an alpha value of
                :math:`\\alpha = \\alpha_0 \\pi`
            #.  The endpoints of the Kaiser-Bessel function tend to
                zero faster than :math:`(1+2\\alpha) / \\exp(\\alpha)`
        Warnings:
            #.  None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are :math:`1/I_0(\\alpha)`.
                For small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha values beyond
                :math:`2.5\\pi` this effect is negligible.
    """
    fname = "diffrec.window_functions.kbal"
    error_check.check_is_real(x, "x", fname)
    error_check.check_positive(W, "W", fname)

    # Window functions have an odd number of points.
    x = 1.0 - np.square(2.0*x/W)

    # Compute window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in
    alpha *= np.pi
    w_func = iv(0.0, alpha * np.sqrt(x))/iv(0.0, alpha)
    return w_func

def kbmdal(x, W, al):
    """
        Purpose:
            Create a modifed Kaiser-Bessel window
            with a user specified alpha parameter.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
            :al (*float*):
                The alpha parameter :math:`\\alpha_0`.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel alpha window of width
                and w_in spacing dx between points.
        Notes:
            #.  The Kaiser-Bessel window is computed using the
                modified Bessel Function of the First Kind. It's
                value is :math:`y = I_0(\\alpha\\sqrt{1-4x^2/w^2})/I_0(\\alpha)`,
                where w is the window width.
            #.  We automatically multiply the alpha parameter by pi,
                so the kbal window function has an alpha value of
                :math:`\\alpha = \\alpha_0\\pi`
            #.  The endpoints of the Kaiser-Bessel function tend to
                zero faster than :math:`(1+2\\alpha)) / \\exp(\\alpha)`
        Warnings:
            #.  None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are :math:`1/I_0(\\alpha)`.
                For small values of alpha this can create Gibb's like
                effects in reconstruction due to the large
                discontinuity in the window. For alpha values beyond
                :math:`2.5\\pi` this effect is negligible.
    """
    fname = "diffrec.window_functions.kbal"
    error_check.check_is_real(x, "x", fname)
    error_check.check_positive(W, "W", fname)

    # Window functions have an odd number of points.
    x = 1.0 - np.square(2.0*x/W)

    # Compute window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in
    alpha *= np.pi
    w_func = (iv(0.0, alpha * np.sqrt(x))-1.0)/(iv(0.0, alpha)-1.0)
    return w_func

# TODO
def window_width(res, normeq, fsky, fres, rho_dot,
                 sigma, bfac=True, Return_P=False):
    """
    Purpose:
        Compute the window width as a function of ring radius.
        This is given from MTR86 Equations 19, 32, and 33.
    Variables:
        :res (*float*):
            The requested resolution.
        :normeq (*float*):
            The normalized equivalent width. Unitless.
        :fsky (*float* or *np.ndarray*):
            The sky frequency.
        :fres (*float* or *np.ndarray*):
            The Fresnel scale.
        :rdot (*float*) or (*np.ndarray*):
            The time derivative of the ring radius.
    Output:
        :w_vals (*np.ndarray*):
            The window width as a function of ring radius.
    """
    if bfac:
        omega = TWO_PI * fsky
        alpha = np.square(omega * sigma) / (2.0 * rho_dot)
        P = res / (alpha * np.square(fres))

        # Create a variable specifying where P>1 occurs.
        Prange = (P > 1.0).nonzero()[0]

        if (np.size(Prange) == 0):
            raise IndexError(
                """
                \r\tError Encountered:
                \r\t\trss_ringoccs.diffrec.DiffractionCorrection\n
                \r\tThe P parameter in window width computation
                \r\tis less than one for the entirety of the
                \r\tdata set. Either rho_dot_km_vals is too small,
                \r\tor F_km_vals is too large. Request a coarser
                \r\tresolution, or check your data for errors.\n
                \r\tAlternatively, you may set bfac=False as
                \r\ta keyword to skip the use of the P parameter.
                \r\tThis may result in inaccurate window width.
                """
            )
        else:
            pass

        P = P[Prange]
        alpha = alpha[Prange]

        w_vals = np.zeros(np.size(rho_dot))
        w_vals[Prange] = special_functions.resolution_inverse(P)/alpha

    else:
        w_vals = 2.0*np.square(fres)/res
        Prange = (fres > 0.0).nonzero()[0]

    w_vals *= normeq

    if Return_P:
        return w_vals, Prange
    else:
        return w_vals

def window_norm(ker, dx, f_scale):
    """
    Purpose:
        Compute the window normalization
    Arguments:
        :ker (*np.ndarray*):
            The Fresnel Kernel.
        :dx (*float*):
            The spacing between points in the window.
            This is equivalent to the sample spacing.
            This value is in kilometers.
        :f_scale (*np.ndarray*):
            The Fresnel Scale in kilometers.
    Outputs:
        :norm_fact (*float*):
            The normalization of the input
            Fresnel Kernel.
    """
    try:
        return _window_functions.window_norm(ker, dx, f_scale)
    except KeyboardInterrupt:
        raise
    except:
        raise TypeError(
            """
            \r\tError: rss_ringoccs
            \r\t\tdiffrec.window_functions.window_norm\n
            \r\tInput should be a numpy array of real numbers (ints or floats),
            \r\tand two positive floating point numbers.
            \r\tUsage:
            \r\t\t>>> x = numpy.arange(-5, 5, 0.1)
            \r\t\t>>> dx = 0.25
            \r\t\t>>> f_scale = 4.0
            \r\t\t>>> y = window_norm(x, dx, f_scale)
            """
        )

func_dict = {
    "rect":   {"func": rect,   "normeq": 1.00000000, "wnum": 0},
    "coss":   {"func": coss,   "normeq": 1.50000000, "wnum": 1},
    "kb20":   {"func": kb20,   "normeq": 1.49634231, "wnum": 2},
    "kb25":   {"func": kb25,   "normeq": 1.65191895, "wnum": 3},
    "kb35":   {"func": kb35,   "normeq": 1.92844639, "wnum": 4},
    "kbmd20": {"func": kbmd20, "normeq": 1.52048382, "wnum": 5},
    "kbmd25": {"func": kbmd25, "normeq": 1.65994438, "wnum": 6},
    "kbmd35": {"func": kbmd35, "normeq": 1.52048382, "wnum": 7}
}
