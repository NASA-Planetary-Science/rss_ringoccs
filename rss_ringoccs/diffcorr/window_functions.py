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
        Function:
            rect
        Purpose:
            Create the rectangular window function
        Variables:
            W:      Window width.
            dx:     Width of one point (Or one bin).
        Outputs:
            w_func: The rectungular window function of width w_in
                    and spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
        Notes:
            The rectangular window function is the unit function.
            That is, it is equal to one across it's entire domain. It
            is used in the Fresnel Inversion and Fresnel Transform
            process to act as a 'hat' function, zeroing out a
            function outside of the domain of rect's definition, and
            having no effect within that domain.
        References:
            [1] https://en.wikipedia.org/wiki/Window_function
            [2] https://en.wikipedia.org/wiki/Rectangular_function
            [3] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        Examples:
            Create a rectangular window function of width 10,
            with 0.01 spacing spacing between points.
                In [1]: import diffcorr as dc
                In [2]: dc.rect(10,0.1)
        History:
            Translated from IDL: RJM - 2018/05/15 9:03 A.M.
            Lowercase variables: RJM - 2018/05/16 1:29 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
        Function:
            coss
        Purpose:
            Create cosine squared window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
        Outputs:
            w_func: The squared cosine window function of width w_in
                    and spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
        History:
            Translated from IDL: RJM - 2018/05/15 9:41 A.M.
            Lowercase variables: RJM - 2018/05/16 1:34 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
        Function:
            kb25
        Purpose:
            Create Kaiser-Bessel 2.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
        Outputs:
            w_func: The Kaiser-Bessel 2.5 window of width w_in and
                    spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
            [3] On approximating the Modified Bessel Function of the
                First Kind and Toader-Qi Mean, Yang, ZH. & Chu, YM.
                J Inequal Appl (2016): 40., Springer,
                https://doi.org/10.1186/s13660-016-0988-1
        History:
            Translated from IDL: RJM - 2018/05/15 9:43 A.M.
            Lowercase variables: RJM - 2018/05/16 3:23 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
        Function:
            kb35
        Purpose:
            Create Kaiser-Bessel 3.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
        Outputs:
            w_func: The Kaiser-Bessel 3.5 window of width w_in and
                    spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
            [3] On approximating the Modified Bessel Function of the
                First Kind and Toader-Qi Mean, Yang, ZH. & Chu, YM.
                J Inequal Appl (2016): 40., Springer,
                https://doi.org/10.1186/s13660-016-0988-1
        History:
            Translated from IDL: RJM - 2018/05/15 9:43 A.M.
            Lowercase variables: RJM - 2018/06/16 3:26 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
        Function:
            kbmd
        Purpose:
            Create the Modified Kaiser-Bessel 2.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
        Outputs:
            w_func: The Modified Kaiser-Bessel 2.5 window of width
                    w_in and spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] https://en.wikipedia.org/wiki/Window_function
        History:
            Translated from IDL: RJM - 2018/05/15 9:44 A.M.
            Lowercase variables: RJM - 2018/05/16 3:34 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
    w_func = (iv(0.0, TWO_PI * np.sqrt((1.0 - 4.0*x*x))) - 1.0)/(IV0_20 - 1.0)
    return w_func

def kbmd25(w_in, dx, error_check=True):
    """
        Function:
            kbmd
        Purpose:
            Create the Modified Kaiser-Bessel 2.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
        Outputs:
            w_func: The Modified Kaiser-Bessel 2.5 window of width
                    w_in and spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] https://en.wikipedia.org/wiki/Window_function
        History:
            Translated from IDL: RJM - 2018/05/15 9:44 A.M.
            Lowercase variables: RJM - 2018/05/16 3:34 P.M.
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
    w_func = (iv(0.0, alpha * np.sqrt((1.0 - 4.0*x*x))) - 1.0)/(IV0_25-1.0)
    return w_func

def kbal(w_in, dx, al, error_check=True):
    """
        Function:
            kb35
        Purpose:
            Create Kaiser-Bessel 3.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
            al:     The alpha parameter.
        Outputs:
            w_func: The Kaiser-Bessel alpha window of width w_in and
                    spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
            [3] On approximating the Modified Bessel Function of the
                First Kind and Toader-Qi Mean, Yang, ZH. & Chu, YM.
                J Inequal Appl (2016): 40., Springer,
                https://doi.org/10.1186/s13660-016-0988-1
        History:
            Created: RJM - 2018/05/16 3:50 P.M.

    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\twindow_functions.kb25:\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\twindow_functions.kb25:\n"
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
        Function:
            kb35
        Purpose:
            Create Kaiser-Bessel 3.5 window.
        Variables:
            W:      Window width.
            dx:     Width of one point.
            al:     The alpha parameter.
        Outputs:
            w_func: The Kaiser-Bessel alpha window of width w_in and
                    spacing dx between points.
        Dependencies:
            [1] diffcorr
            [2] sys
            [3] numpy
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
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
            [2] https://en.wikipedia.org/wiki/Window_function
            [3] On approximating the Modified Bessel Function of the
                First Kind and Toader-Qi Mean, Yang, ZH. & Chu, YM.
                J Inequal Appl (2016): 40., Springer,
                https://doi.org/10.1186/s13660-016-0988-1
        History:
            Created: RJM - 2018/05/16 3:50 P.M.

    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    """
    if error_check:
        if not isinstance(w_in, float):
            try:
                w_in = float(w_in)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\twindow_functions.kb25:\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
                "\twindow_functions.kb25:\n"
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
        Function:
            window_width
        Purpose:
            Compute the window width as a function of ring radius.
            This is given from MTR86 Equations 19, 32, and 33.
        Variables:
            res:
                The requested resolution.
            normeq:
                The normalized equivalent width. Unitless.
            fsky:
                The sky frequency.
            fres:
                The Fresnel scale.
            rdot:
                The time derivative of the ring radius.
        Output:
            w_vals:  The window width as a function of ring radius.
        History:
            Translated from IDL: RJM - 2018/05/15 8:38 P.M.
    """
    if error_check:
        if (not isinstance(res, float)):
            try:
                res = float(res)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tFirst input must be positive.\n"
                    "\t\tYour input : %f\n"
                    % (res)
            )
        else:
            pass

        if (not isinstance(normeq, float)):
            try:
                res = float(ref)
            except (TypeError, ValueError):
                raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfsky must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            elif (np.min(fsky) < 0.0):
                raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\tfres must be a positive float or\n"
                    "\t\ta numpy array of positive floats.\n"
                )
            elif (np.min(fres) < 0.0):
                raise ValueError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
                    "\twindow_functions.window_width:\n"
                    "\t\trdot must be a float or a numpy\n"
                    "\t\tarray of floats.\n"
                )
            else:
                pass

        except (ValueError, TypeError) as errmes:
            raise TypeError(
                    "\n\tError Encountered:\n"
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                    "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
                "\trss_ringoccs: Diffcorr Subpackage\n"
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
    # Freespace Integral
    T1 = np.abs(np.sum(ker) * dx)

    # Normalization Factor
    norm_fact = SQRT_2 * f_scale / T1
    return norm_fact

def get_range_actual(rho, rng, w_vals):
    """
        Function:
            get_range_actual
        Purpose:
            Compute the possible allowed range for processing, taking
            into consideration available data (rho) and the requested region.
        Variables:
            RHO:        Radial range of the data.
            RANGE:      Requested start/end points for processing.
            W_KM_VALS:  Window width as a function of ring radius.
        Output:
            START:  The allowed starting point for processing.
            N_USED: The number of points allowed for processing.
        History:
            Translated from IDL: RJM - 2018/05/15 3:19 P.M.
    """
    if (not check_real(rho)):
        sys.exit("Rho must be an array of real numbers")
    if (np.min(rho) < 0.0):
        sys.exit("Rho must be positive")
    if (np.size(rng) != 2):
        sys.exit("Range must have format rng = [a,b]")
    if (not check_pos_real(np.min(rng))):
        sys.exit("Range must be positive")
    if (not check_real(w_vals)):
        sys.exit("w_vals must be real")
    if (np.min(w_vals) < 0.0): 
        sys.exit("w_vals must be positive")
    if (np.min(rng) > np.max(rho)):
        sys.exit("Requested range GREATER than available data.")
    if (np.max(rng) < np.min(rho)):
        sys.exit("Requested range LESS than available data.")
    w_max       = np.max(w_vals)
    rho_min_lim = np.min(rho)+np.ceil(w_max/2.0)
    rho_max_lim = np.max(rho)-np.ceil(w_max/2.0)
    rho_start   = rho[np.min((rho >= np.min(rng)).nonzero())]
    rho_end     = rho[np.max((rho <= np.max(rng)).nonzero())]
    rho_min     = np.max([rho_min_lim,rho_start])
    rho_max     = np.min([rho_max_lim,rho_end])
    start       = int(np.min((rho >= rho_min).nonzero()))
    finish      = int(np.max((rho <= rho_max).nonzero()))
    n_used      = 1 + (finish - start)
    return start, n_used

func_dict = {
    "rect":   {"func": rect,   "normeq": 1.00000000},
    "coss":   {"func": coss,   "normeq": 1.50000000},
    "kb20":   {"func": kb20,   "normeq": 1.49634231},
    "kb25":   {"func": kb25,   "normeq": 1.65191895},
    "kb35":   {"func": kb35,   "normeq": 1.92844639},
    "kbmd20": {"func": kbmd20, "normeq": 1.52048174},
    "kbmd25": {"func": kbmd25, "normeq": 1.65994218}
}