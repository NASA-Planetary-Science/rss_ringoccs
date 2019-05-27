"""
    Purpose:
        Provide a suite of window functions and
        functions related to the normalized equivalent
        width of a given array.
    Dependencies:
        #. numpy
        #. spicy
"""

import numpy as np
from scipy.special import lambertw, iv
try:
    from rss_ringoccs._ufuncs import _window_functions
except (ImportError, ModuleNotFoundError):
    print(
        """
            Error: rss_ringoccs.diffrec.window_functions
            \tCould Not Import C Code. Stricly Using Python Code.
            \tThis is signicantly slower. There was most likely an error
            \tin your installation of rss_ringoccs. To use the C Code,
            \tdownload a C Compiler (GCC) and see the User's Guide for
            \tinstallation instructions.
        """
    )

# Declare constants for multiples of pi.
TWO_PI = 6.283185307179586476925287
SQRT_PI_2 = 1.253314137315500251207883
RCP_SQRT_2 = 0.7071067811865476

# Declare constants for the reciprocal of e and the square root of 2.
RCPR_E = 0.3678794411714423215955238
SQRT_2 = 1.414213562373095048801689

# Declare constants for various Bessel function inputs (I_0(x)).
IV0_20 = 87.10850209627940
IV0_25 = 373.02058499037486
IV0_35 = 7257.7994923041760

def rect(x, W):
    """
        Purpose:
            Compute the rectangular window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.rect(x, W)
    except (TypeError, ValueError, NameError):
        w_func = np.zeros(np.size(x)) + 1.0
        w_func[(np.abs(x) >= W/2.0).nonzero()] = 0.0
        return w_func

def coss(x, W):
    """
        Purpose:
            Compute the squared cosine window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.coss(x, W)
    except (TypeError, ValueError, NameError):
        x *= np.pi/W

        w_func = np.square(np.cos(x))
        w_func[(np.abs(x) >= W/2.0).nonzero()] = 0.0
        return w_func

def kb20(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 2.0 window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kb20(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Compute window function.
        w_func = (iv(0.0, TWO_PI*np.sqrt(x)))/IV0_20
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kb25(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 2.5 window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kb25(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Alpha value for kb25 is 2.5.
        alpha = 2.5 * ONE_PI

        # Compute window function.
        w_func = (iv(0.0, alpha*np.sqrt(x)))/IV0_25
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kb35(x, W):
    """
        Purpose:
            Compute the Kaiser-Bessel 3.5 window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kb35(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Alpha value for kb35 is 3.5.
        alpha = 3.5 * ONE_PI

        # Compute window function.
        w_func = (iv(0.0, alpha*np.sqrt(x)))/IV0_35
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kbmd20(x, W):
    """
        Purpose:
            Compute the Modifed Kaiser-Bessel 2.0
            window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kbmd20(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Compute window function.
        w_func = (iv(0.0, TWO_PI*np.sqrt(x)) - 1.0)/(IV0_20 - 1.0)
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kbmd25(x, W):
    """
        Purpose:
            Compute the Modifed Kaiser-Bessel 2.5
            window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kbmd25(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Alpha value for kbmd25 is 2.5.
        alpha = 2.5 * ONE_PI

        # Compute window function.
        w_func = (iv(0.0, alpha*np.sqrt(x)) - 1.0)/(IV0_25 - 1.0)
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kbmd35(x, W):
    """
        Purpose:
            Compute the Modifed Kaiser-Bessel 3.5
            window function.
        Arguments:
            :w_in (*float*):
                The width of the window function to be computed.
                This value is in kilometers.
            :dx (*float*):
                The spacing between points in the window.
                This is equivalent to the sample spacing.
                This value is in kilometers.
        Outputs:
            :w_func (*np.ndarray*):
                Window function of width
                w_in and with sample spacing dx.
    """
    try:
        return _window_functions.kbmd35(x, W)
    except (TypeError, ValueError, NameError):
        x = 1.0 - np.square(2.0*x/W)

        # Alpha value for kbmd35 is 3.5.
        alpha = 3.5 * ONE_PI

        # Compute window function.
        w_func = (iv(0.0, alpha*np.sqrt(x)) - 1.0)/(IV0_35 - 1.0)
        w_func[(x <= 0).nonzero()] = 0.0
        return w_func

def kbal(w_in, dx, al, error_check=True):
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
                The alpha parameter.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel alpha window of width
                w_in and spacing dx between points.
        Notes:
            [1] The Kaiser-Bessel window is computed using the 
                modified Bessel Function of the First Kind. It's
                value is y = I_0(alpha*sqrt(1-4x^2/w^2))/I_0(alpha),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kbal window function has an alpha value of
                alpha = al * pi
            [3] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1+2 * alpha)) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1/I_0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha values beyond
                2.5 * pi this effect is negligible.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbal:\n"
                    "\t\tFirst input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(w_in).__name__)
                )
        else:
            pass
    
        if (w_in <= 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbal:\n"
                "\t\tFirst input must positive.\n"
                "\t\tYour input: %f\n"
                % (w_in)
            )
        else:
            pass

        if not isinstance(dx, float):
            try:
                dx = float(dx)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbal:\n"
                    "\t\tSecond input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(dx).__name__)
                )
        else:
            pass

        if (dx <= 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbal:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

        if not isinstance(al, float):
            try:
                al = float(al)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbal:\n"
                    "\t\tThird input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(al).__name__)
                )
        else:
            pass

        if (al < 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbal:\n"
                "\t\tSecond input must non-negative.\n"
                "\t\tYour input: %f\n"
                % (al)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in
    alpha  =  al * np.pi
    w_func = iv(0.0, alpha * np.sqrt((1.0 - 4.0*x*x)))/iv(0.0, alpha)
    return w_func

def kbmdal(w_in, dx, al, error_check=True):
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
                The alpha parameter.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel alpha window of width
                and w_in spacing dx between points.
        Notes:
            [1] The Kaiser-Bessel window is computed using the 
                modified Bessel Function of the First Kind. It's
                value is y = I_0(alpha*sqrt(1-4x^2/w^2))/I_0(alpha),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kbal window function has an alpha value of
                alpha = al * pi
            [3] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1+2 * alpha)) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1/I_0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha values beyond
                2.5 * pi this effect is negligible.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbmdal:\n"
                    "\t\tFirst input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(w_in).__name__)
                )
        else:
            pass
    
        if (w_in <= 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbmdal:\n"
                "\t\tFirst input must positive.\n"
                "\t\tYour input: %f\n"
                % (w_in)
            )
        else:
            pass

        if not isinstance(dx, float):
            try:
                dx = float(dx)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbmdal:\n"
                    "\t\tSecond input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(dx).__name__)
                )
        else:
            pass

        if (dx <= 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbmdal:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

        if not isinstance(al, float):
            try:
                al = float(al)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbmdal:\n"
                    "\t\tThird input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(al).__name__)
                )
        else:
            pass

        if (al < 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.kbmdal:\n"
                "\t\tSecond input must non-negative.\n"
                "\t\tYour input: %f\n"
                % (al)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in
    alpha  =  al * np.pi
    w_func = (iv(0.0, alpha*np.sqrt((1.0-4.0*x*x)))-1.0)/(iv(0.0, alpha)-1.0)
    return w_func

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
                "\n\tError Encountered:\n"
                "\t\trss_ringoccs.diffrec.DiffractionCorrection\n\n"
                "\tThe P parameter in window width computation\n"
                "\tis less than one for the entirety of the\n"
                "\tdata set. Either rho_dot_km_vals is too small,\n"
                "\tor F_km_vals is too large. Request a coarser\n"
                "\tresolution, or check your data for errors.\n\n"
                "\tAlternatively, you may set bfac=False as\n"
                "\ta keyword to skip the use of the P parameter.\n"
                "\tThis may result in inaccurate window widths."
            )
        else:
            pass

        P = P[Prange]
        alpha = alpha[Prange]
        P1 = P/(1-P)
        P2 = P1*np.exp(P1)

        # LambertW returns nans far values close to zero, so round this.
        P2 = np.around(P2, decimals=16)

        # For values near -1/e, LambertW(x) is roughly -1. 
        crange1 = ((RCPR_E + P2) < 1.0e-16).nonzero()[0]
        crange2 = ((RCPR_E + P2) >= 1.0e-16).nonzero()[0]
        w_vals = np.zeros(np.size(rho_dot))

        if (np.size(crange1) > 0):
            w_vals[Prange[crange1]] = 2.0*np.square(fres[Prange[crange1]])/res
        else:
            pass

        if (np.size(crange2) > 0):
            w_vals[Prange[crange2]] = np.abs(lambertw(P2)-P1)/alpha
        else:
            pass
    else:
        w_vals = 2.0*np.square(fres)/res

    w_vals *= normeq

    if Return_P:
        return w_vals, Prange
    else:
        return w_vals

def normalize(dx, ker, f_scale):
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
    # Freespace Integral
    T1 = np.abs(np.sum(ker) * dx)

    # Normalization Factor
    norm_fact = SQRT_2 * f_scale / T1
    return norm_fact

func_dict = {
    "rect":   {"func": rect,   "normeq": 1.00000000},
    "coss":   {"func": coss,   "normeq": 1.50000000},
    "kb20":   {"func": kb20,   "normeq": 1.49634231},
    "kb25":   {"func": kb25,   "normeq": 1.65191895},
    "kb35":   {"func": kb35,   "normeq": 1.92844639},
    "kbmd20": {"func": kbmd20, "normeq": 1.52048382},
    "kbmd25": {"func": kbmd25, "normeq": 1.65994438},
    "kbmd35": {"func": kbmd35, "normeq": 1.52048382}
}