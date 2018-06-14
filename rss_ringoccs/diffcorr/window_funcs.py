from .misc_tools import *
import sys
import numpy as np
from scipy.special import iv

def rect(w_in, dx):
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
    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    w_func = np.zeros(nw_pts) + 1.0
    return w_func

def coss(w_in, dx):
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
    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
    w_func = np.cos(np.pi * x / w_in)**2
    return w_func

def kb25(w_in, dx):
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
    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
    alpha  = 2.5*np.pi
    w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
    return w_func

def kb35(w_in, dx):
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
    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
    alpha  = 3.5 * np.pi
    w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
    return w_func

def kbal(w_in, dx, al):
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
    """
    tw  = check_pos_real(w_in)
    tdx = check_pos_real(dx)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
    alpha  =  al * np.pi
    w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
    return w_func

def kbmd(w_in, dx):
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
    tdx = check_pos_real(dx)
    tw  = check_pos_real(w_in)
    if (not tdx) or (not tw):
        sys.exit("Input must be two positive real numbers")
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
    alpha  = 2.5 * np.pi
    w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
    return w_func

func_dict = {
    "rect" : {"func" : rect, "normeq" : 1.00000000},
    "coss" : {"func" : coss, "normeq" : 1.50000000},
    "kb25" : {"func" : kb25, "normeq" : 1.65191895},
    "kb35" : {"func" : kb35, "normeq" : 1.92844639},
    "kbmd" : {"func" : kbmd, "normeq" : 1.65994218}
    }