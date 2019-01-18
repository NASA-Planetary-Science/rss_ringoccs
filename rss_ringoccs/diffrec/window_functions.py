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

# Declare constants for multiples of pi.
TWO_PI = 6.283185307179586476925287
ONE_PI = 3.141592653589793238462643

# Declare constants for the reciprocal of e and the square root of 2.
RCPR_E = 0.3678794411714423215955238
SQRT_2 = 1.414213562373095048801689

# Declare constants for various Bessel function inputs (I_0(x)).
IV0_20 = 87.10850209627940
IV0_25 = 373.02058499037486
IV0_35 = 7257.7994923041760

def rect(w_in, dx, error_check=True):
    """
        Purpose:
            Create the rectangular window function
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point (Or one bin).
        Outputs:
            :w_func (*np.ndarray*):
                The rectungular window function of width w_in
                and spacing dx between points. Numpy array.
        Notes:
            The rectangular window function is the unit function.
            That is, it is equal to one across it's entire domain. It
            is used in the Fresnel Inversion and Fresnel Transform
            process to act as a 'hat' function, zeroing out a
            function outside of the domain of rect's definition, and
            having no effect within that domain.
        Examples:
            Create a rectangular window function of width 10,
            with 0.01 spacing spacing between points.
                In [1]: from rss_ringoccs import diffrec
                In [2]: diffrec.rect(10, 0.1)
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.rect:\n"
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
                "\twindow_functions.rect:\n"
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
                    "\twindow_functions.rect:\n"
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
                "\twindow_functions.rect:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    w_func = np.zeros(nw_pts) + 1.0
    return w_func

def coss(w_in, dx, error_check=True):
    """
        Purpose:
            Create cosine squared window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            :w_func (*np.ndarray*):
                The squared cosine window function of width w_in
                and spacing dx between points.
        Notes:
            [1] This window function is defined as
                y = cos^2(x * pi/w), where w is the window width.
            [2] The squared cosine window has the advantage of
                evaluating to zero at the endpoints of the window,
                meaning Gibb's effects from discontinuities do not
                arise like in the rectangular window and, to a much
                lesser degree, the Kaiser-Bessel windows.
            [3] The squared cosine window is wider than the all of
                the Kaiser-Bessel windows for alpha > 2.5*pi. Alpha
                is a paramter in the Kaiser-Bessel window.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.coss:\n"
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
                "\twindow_functions.coss:\n"
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
                    "\twindow_functions.coss:\n"
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
                "\twindow_functions.coss:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # The number of points in a window must be odd.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Set the independent variable for the window.
    x  = ONE_PI * (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    w_func = np.cos(x)*np.cos(x)
    return w_func

def kb20(w_in, dx, error_check=True):
    """
        Purpose:
            Create Kaiser-Bessel 2.0 window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            :w_func (*np.ndarray*):
                The squared cosine window function of width w_in
                and spacing dx between points.
        Notes:
            [1] This window function is defined as
                y = cos^2(x * pi/w), where w is the window width.
            [2] The squared cosine window has the advantage of
                evaluating to zero at the endpoints of the window,
                meaning Gibb's effects from discontinuities do not
                arise like in the rectangular window and, to a much
                lesser degree, the Kaiser-Bessel windows.
            [3] The squared cosine window is wider than the all of
                the Kaiser-Bessel windows for alpha > 2.5*pi. Alpha
                is a paramter in the Kaiser-Bessel window.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kb20:\n"
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
                "\twindow_functions.kb20:\n"
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
                    "\twindow_functions.kb20:\n"
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
                "\twindow_functions.kb20:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    # Compute window function.
    w_func = iv(0.0, TWO_PI * np.sqrt((1.0 - 4.0*x*x))) / IV0_20
    return w_func

def kb25(w_in, dx, error_check=True):
    """
        Purpose:
            Create Kaiser-Bessel 2.5 window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel 2.5 window of width w_in
                and spacing dx between points.
        Notes:
            [1] The Kaiser-Bessel window is computed using the 
                modified Bessel Function of the First Kind. It's
                value is y = I_0(alpha*sqrt(1-4x^2/w^2))/I_0(alpha),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kb25 window function has an alpha value of
                alpha = 2.5 * pi
            [3] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1+2 * alpha)) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1/I_0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha = 2.5 * pi,
                the endpoint is 0.00268082, which is close to zero.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kb25:\n"
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
                "\twindow_functions.kb25:\n"
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
                    "\twindow_functions.kb25:\n"
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
                "\twindow_functions.kb25:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    # Alpha value for kb25 is 2.5.
    alpha = 2.5 * ONE_PI

    # Compute window function.
    w_func = iv(0.0, alpha * np.sqrt((1.0 - 4.0*x*x))) / IV0_25
    return w_func

def kb35(w_in, dx, error_check=True):
    """
        Purpose:
            Create Kaiser-Bessel 3.5 window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            :w_func (*np.ndarray*):
                The Kaiser-Bessel 3.5 window of width w_in
                and spacing dx between points.
        Notes:
            [1] The Kaiser-Bessel window is computed using the 
                modified Bessel Function of the First Kind. It's
                value is y = I_0(alpha*sqrt(1-4x^2/w^2))/I_0(alpha),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kb35 window function has an alpha value of
                alpha = 3.5 * pi
            [3] The endpoints of the Kaiser-Bessel function tend to
                zero faster than (1+2 * alpha)) / exp(alpha)
        Warnings:
            [1] None of the Kaiser-Bessel windows evaluate to zero at
                their endpoints. The endpoints are 1/I_0(alpha). For
                small values of alpha this can create Gibb's like
                effects in reconstruction do to the large
                discontinuity in the window. For alpha = 3.5 * pi,
                the endpoint is 0.000137783, which is close to zero.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kb35:\n"
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
                "\twindow_functions.kb35:\n"
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
                    "\twindow_functions.kb35:\n"
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
                "\twindow_functions.kb35:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    # Alpha value for kb35 is 3.5.
    alpha = 3.5 * ONE_PI

    # Compute window function.
    w_func = iv(0.0, alpha * np.sqrt((1.0 - 4.0*x*x))) / IV0_35
    return w_func

def kbmd20(w_in, dx, error_check=True):
    """
        Purpose:
            Create the Modified Kaiser-Bessel 2.5 window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            w_func (*np.ndarray*):
                The Modified Kaiser-Bessel 2.5 window of
                width w_in and spacing dx between points.
        Notes:
            [1] The Modified Kaiser-Bessel window is computed using
                the modified Bessel Function of the First Kind. It's
                value is:
                y = (I_0(alpha*sqrt(1-4x^2/w^2))-1)/(I_0(alpha)-1),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kbmd window function has an alpha value of
                alpha = 2.5 * pi
            [3] The endpoints of the Modified Kaiser-Bessel function
                are zero. That means that the modified version has
                no discontinuities in it.
            [4] The Kaiser-Bessel functions and the modified
                Kaiser-Bessel functions are equal at the center of
                the window.
        Warnings:
            [1] Unlike the Kaiser-Bessel functions, the modified
                Kaiser-Bessel functions evaluate to zero at the
                endpoints of the window.
            [2] Small alpha values will result in the Kaiser-Bessel
                function and the modified Kaiser-Bessel disagreeing
                dramatically. For example, alpha = 0 gives a constant
                curve for the Kaiser-Bessel, but a bell-shaped curve
                for the modified version.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbmd20:\n"
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
                "\twindow_functions.kbmd20:\n"
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
                    "\twindow_functions.kbmd20:\n"
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
                "\twindow_functions.kbmd20:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    # Compute window function.
    w_func = (iv(0.0, TWO_PI * np.sqrt((1.0 - 4.0*x*x))) - 1.0)/(IV0_20 - 1.0)
    return w_func

def kbmd25(w_in, dx, error_check=True):
    """
        Purpose:
            Create the Modified Kaiser-Bessel 2.5 window.
        Variables:
            :W (*float*):
                Window width.
            :dx (*float*):
                Width of one point.
        Outputs:
            w_func (*np.ndarray*):
                The Modified Kaiser-Bessel 2.5 window of
                width w_in and spacing dx between points.
        Notes:
            [1] The Modified Kaiser-Bessel window is computed using
                the modified Bessel Function of the First Kind. It's
                value is:
                y = (I_0(alpha*sqrt(1-4x^2/w^2))-1)/(I_0(alpha)-1),
                where w is the window width.
            [2] We automatically multiply the alpha parameter by pi,
                so the kbmd window function has an alpha value of
                alpha = 2.5 * pi
            [3] The endpoints of the Modified Kaiser-Bessel function
                are zero. That means that the modified version has
                no discontinuities in it.
            [4] The Kaiser-Bessel functions and the modified
                Kaiser-Bessel functions are equal at the center of
                the window.
        Warnings:
            [1] Unlike the Kaiser-Bessel functions, the modified
                Kaiser-Bessel functions evaluate to zero at the
                endpoints of the window.
            [2] Small alpha values will result in the Kaiser-Bessel
                function and the modified Kaiser-Bessel disagreeing
                dramatically. For example, alpha = 0 gives a constant
                curve for the Kaiser-Bessel, but a bell-shaped curve
                for the modified version.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.kbmd25:\n"
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
                "\twindow_functions.kbmd25:\n"
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
                    "\twindow_functions.kbmd25:\n"
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
                "\twindow_functions.kbmd25:\n"
                "\t\tSecond input must positive.\n"
                "\t\tYour input: %f\n"
                % (dx)
            )
        else:
            pass

    # Window functions have an odd number of points.
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)

    # Compute argument of window function.
    x = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx / w_in

    # Alpha value for kb25 is 2.5.
    alpha = 2.5 * ONE_PI

    # Compute window function.
    w_func = (iv(0.0, alpha * np.sqrt((1.0 - 4.0*x*x))) - 1.0)/(IV0_25-1.0)
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
    alpha  =  al * ONE_PI
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
    alpha  =  al * ONE_PI
    w_func = (iv(0.0, alpha*np.sqrt((1.0-4.0*x*x)))-1.0)/(iv(0.0, alpha)-1.0)
    return w_func

def window_width(res, normeq, fsky, fres, rdot,
                 sigma = 2.e-13, bfac=True, error_check=True):
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
    if error_check:
        if (not isinstance(res, float)):
            try:
                res = float(res)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tFirst input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(res).__name__)
                )
        else:
            pass
        
        if (res <= 0):
            raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tFirst input must be positive.\n"
                    "\t\tYour input : %f\n"
                    % (res)
            )
        else:
            pass

        if (not isinstance(normeq, float)):
            try:
                res = float(res)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tSecond input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(normeq).__name__)
                )
        else:
            pass
        
        if (normeq <= 0):
            raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tSecond input must be positive.\n"
                    "\t\tYour input : %f\n"
                    % (normeq)
            )
        else:
            pass
        
        try:
            fsky = np.array(fsky)
            if not (np.all(np.isreal(fsky))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfsky must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            elif (np.min(fsky) < 0.0):
                raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfsky must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            else:
                pass

            fres = np.array(fres)
            if not (np.all(np.isreal(fres))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfres must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            elif (np.min(fres) < 0.0):
                raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfres must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            else:
                pass

            rdot = np.array(rdot)
            if not (np.all(np.isreal(rdot))):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\trdot must be a float or a numpy\n"
                    "\t\tarray of floats.\n"
                )
            else:
                pass

        except (ValueError, TypeError) as errmes:
            raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tOne of your inputs was unable to\n"
                    "\t\tbe converted into a numpy array.\n"
                    "\tOrginal Error Message: %s\n"
                    % (errmes)
            )
        
        if (not isinstance(sigma, float)):
            try:
                sigma = float(sigma)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tKeyword sigma must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(sigma).__name__)
                )
        else:
            pass

        if (sigma <= 0.0):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.window_width:\n"
                "\t\tKeyword sigma must be positive.\n"
                "\t\tYour input: %f\n"
                % (sigma)
            )
        else:
            pass

        if (not isinstance(bfac, bool)):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.window_width:\n"
                "\t\tbfac keyword must be a boolean.\n"
                "\t\tSet bfac=True or bfac=False.\n"
            )
        else:
            pass

    f2 = fres*fres
    if bfac:
        omega = TWO_PI * fsky
        alpha = omega*omega * sigma*sigma / (2.0 * rdot)
        P = res / (alpha * f2)
        # The inverse exists only if P>1.
        if (np.min(P) <= 1.0):
            raise ValueError(
                "\n\tWarning: Bad Points!\n"
                "\tEither rho_dot_kms_vals, F_km_vals, or res_km is\n"
                "\ttoo small, or sigma is too large. It is also\n"
                "\tpossible that the frequency is too large.\n"
                "\tRequest coarser resolution, check your input\n"
                "\tdata, or set bfac=False as a keyword.\n"
            )
        else:
            pass

        P1 = P/(1-P)
        P2 = P1*np.exp(P1)
        crange1 = ((RCPR_E + P2) < 1.0e-16).nonzero()
        crange2 = ((RCPR_E + P2) >= 1.0e-16).nonzero()
        w_vals = np.zeros(np.size(fres))
        if (np.size(crange1) > 0):
            w_vals[crange1] = 2.0*f2/res
        else:
            pass

        if (np.size(crange2) > 0):
            w_vals[crange2] = np.abs(lambertw(P2)-P1)/alpha
        else:
            pass

        del omega, alpha, P, P1, P2, crange1, crange2
    else:
        w_vals = 2.0*f2/res

    w_vals *= normeq

    return w_vals

def normalize(dx, ker, f_scale, error_check=True):
    if error_check:
        if (not isinstance(dx, float)):
            try:
                dx = float(dx)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.normalize:\n"
                    "\t\tFirst input must be a floating point number.\n"
                    "\t\tYour input has type: %s\n"
                    % (type(dx).__name__)
                )
        else:
            pass

        if (dx <= 0.0):
            raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: diffrec Subpackage\n"
                    "\twindow_functions.normalize:\n"
                    "\t\tFirst input must be positive.\n"
                    "\t\tYour input: %f\n"
                    % (dx)
            )
        else:
            pass
        
        try:
            ker = np.array(ker)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.normalize:\n"
                "\t\tSecond argument could not be converted\n"
                "\t\tinto a numpy array.\n"
            )

        try:
            f_scale = np.array(f_scale)
        except (ValueError, TypeError):
            raise TypeError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.normalize:\n"
                "\t\tThird argument could not be converted\n"
                "\t\tinto a numpy array.\n"
            )
        
        if (not np.all(np.isreal(f_scale))):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.normalize:\n"
                "\t\tThird argument must be real valued.\n"
            )
        elif (np.min(f_scale) <= 0.0):
            raise ValueError(
                "\n\tError Encountered:\n"
                "\trss_ringoccs: diffrec Subpackage\n"
                "\twindow_functions.normalize:\n"
                "\t\tThird argument must be positive.\n"
            )
        else:
            pass
    else:
        pass

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
    "kbmd20": {"func": kbmd20, "normeq": 1.52048174},
    "kbmd25": {"func": kbmd25, "normeq": 1.65994218}
}