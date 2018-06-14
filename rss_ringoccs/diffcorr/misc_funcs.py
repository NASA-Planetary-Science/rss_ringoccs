import numpy as np
from scipy.special import lambertw
from scipy.special import erf
from .misc_tools import *
from .window_funcs import *

region_dict = {
    'maxwell'       : [87410.0,87610.0],
    'maxwellringlet': [87410.0,87610.0],
    'titan'         : [77870.0,77930.0],
    'titanringlet'  : [77870.0,77930.0],
    'huygens'       : [117650.0,117950.0],
    'huygensringlet': [117650.0,117950.0],
    'encke'         : [132900.0,134200.0],
    'enckegap'      : [132900.0,134200.0],
    'all'           : [65000.0,145000.0]
    }

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
    tw = type(w_func)
    if not (tw in vtypes.rtypes):
        sys.exit("Input must be real valued")
    nw      = np.size(w_func)
    normeq  = nw*(np.sum(w_func**2)) / (np.sum(w_func)**2)
    return normeq

def resolution_inverse(x):
    """
        Function:
            resolution_inverse
        Purpose:
            Compute the inverse of y = x/(exp(-x)+x-1)
        Variables:
            x:      A real or complex number.
        Outputs:
            f:      The inverse of x/(exp(-x)+x-1)
        Dependencies:
            [1] numpy
            [2] scipy.special
            [3] sys
        Notes:
            If the input is not in the atypes list, this function
            will execute a system exit, returning the user to
            whichever shell they ran this function from. This
            function wants ints, floats, complex numbers, or arrays.
        Method:
            The inverse of x/(exp(-x)+x-1) is computed using the
            LambertW function. This function is the inverse of
            y = x * exp(x). This is computed using the scipy.special
            lambertw function.
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
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: x = np.array(range(0,1001))*0.001+1.01
                In [4]: y = dc.resolution_inverse(x)
                In [5]: import matplotlib.pyplot as plt
                In [6]: plt.show(plt.plot(x,y))
                (Beautiful plots appear here)
        History:
            Translated from IDL: RJM - 2018/05/15 1:49 P.M.
    """
    if (type(x) not in vtypes.atypes):
        sys.exit("Input must be real or complex")
    f = ((x - 1) * lambertw(np.exp(x / (1 - x)) * x / (1 - x)) + x) / (x - 1)
    return f

def get_range_request(rngreq):
    """
        Function:
            get_range_request
        Purpose:
            This function takes a variety of input ranges and then
            converts them into the form [start,finish].
        Variables:
            rngreq: The start/end points for diffraction correction.
                    Preferred input is rngreq = [a,b]. Arrays are
                    allowed and the range will be set as:
                        rng = [MIN(array),MAX(array)]
                    Finally, certain strings containing a few of the
                    regions of interests within the rings of Saturn
                    are allowed. Permissible strings are:
                        'maxwell', 'titan', 'huygens', and 'encke'.
                    Strings are neither case nor space sensitive.
                    For other planets use rng = [a,b]
        Outputs:
            rng:    A numpy array of the form [start,finish]
        Dependencies:
            [1] numpy
        Notes:
            [1] Inputs should be in kilometers. Numpy arrays, lists,
                and strings are allowed.
            [2] Strings are neither case nor space sensitive.
            [3] If you enter an illegal string, you will be shown a
                list of allowed strings and prompted for one of them.
            [4] If you enter an illegal value you will be prompted
                for a start and end point.
        References:
            [1] A.F. Cook, F.A. Franklin, F.D. Palluconi,
                Saturn's rings—A survey,
                Icarus, Volume 18, Issue 2, 1973, Pages 317-337,
                https://doi.org/10.1016/0019-1035(73)90214-5
            [2] Pollack, J.B., The Rings of Saturn,
                Space Science Reviews, Volume 18, Issue 1, 1975,
                Pages 3–93, https://doi.org/10.1007/BF00350197
        Examples:
            Try some inputs:
                In [1]: import diffcorr as dc
                In [2]: dc.get_range_request('maxwell')
                Out[2]: array([87410., 87610.])
                In [3]: dc.get_range_request('bob')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: 'titan'
                Out[3]: array([77870., 77930.])
            You don't need to wrap your selection in single quotes.
                In [4]: dc.get_range_request('carl')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: huygens
                Out[4]: array([117650., 117950.])
            You can use caps and double quotes too!
                In [5]: dc.get_range_request('jeff')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: "ENCKE"
                Out[5]: array([132900., 134200.])
            Spacing doesn't matter!
                In [6]: dc.get_range_request('george')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: m      aX     w el L
                Out[6]: array([87410., 87610.])
                In [7]: dc.get_range_request([123,456])
                Out[7]: array([123., 456.])
                In [8]: dc.get_range_request(6)
                You only provided one range value. I need two.
                Enter starting radius (in Kilometers): 6
                Enter final radius (in Kilometers): 10
                Out[8]: array([ 6., 10.])
        History:
            Translated from IDL: RJM - 2018/05/15 2:03 P.M.
    """
    ti = type(rngreq)
    if (not ti in vtypes.rtypes) and (ti != type([1])) and (ti != type("Hi")):
        print("I don't understand your requested range.")
        start  = input("Enter starting radius (in Kilometers): ")
        start  = float(start)
        finish = input("Enter final radius (in Kilometers): ")
        finish = float(finish)
        rng    = np.array([start,finish])
    elif (ti == type('Hello')):
        reg = rngreq.replace(" ","").lower()
        if (reg in region_dict):
            rng = np.array(region_dict[reg])
        else:
            print("Illegal range request. Allowed strings are:")
            print("'maxwell', 'titan', 'huygens', 'encke', or 'all'")
            rngreq = input("Please enter requested range string: ")
            rngreq = rngreq.replace(" ","").lower()
            rngreq = rngreq.replace("'","")
            rngreq = rngreq.replace('"',"")
            rng    = np.array(region_dict[rngreq])
    elif (np.size(rngreq) < 2):
        print("You only provided one range value. I need two.")
        start  = input("Enter starting radius (in Kilometers): ")
        start  = float(start)
        finish = input("Enter final radius (in Kilometers): ")
        finish = float(finish)
        rng    = np.array([start,finish])
    else:
        rng = np.array([np.float(np.min(rngreq)),np.float(np.max(rngreq))])
    return rng

def get_norm_eq(wtype):
    """
        Function:
            get_norm_eq
        Purpose:
            Compute the Normalized Equivalent Width from a given WTYPE.
        Variables:
            wtype:      The name of the window used for processing.
        Output:
            norm_eq:    The normalized equivalent width of wtype.
        Dependencies:
            [1] numpy
        Notes:
            [1] The normalized equivalent width is pre-computed using
                a window with infinite resolution. That is, the
                intergral is perform on the functions exactly, not
                approximated with Riemann sums. This value is more
                accurate than the method used in the compute_norm_eq
                function, however only a few select windows are
                included. These windows are:
                    [1] Rectangular Window..................'rect'
                    [2] Squared Cosine Window...............'coss'
                    [3] Kaiser-Bessel 2.5 Window............'kb25'
                    [4] Kaiser-Bessel 3.5 Window............'kb35'
                    [5] Modified Kaiser-Bessel 2.5 Window...'kbmd'
            [2] The numerical values that this function returns are
                slightly different than those quoted in the reference
                below. The MTR86 paper only evaluates to two decimals
                whereas we have done double precision to 8 decimals.
            [3] The input must be a string. The function is neither
                case nor space sensitive.
            [4] The normalized equivalent width is a unitless value.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        Examples:
            Test some input strings.
                In [1]: import diffcorr as dc
                In [2]: dc.get_norm_eq('kb25')
                Out[2]: 1.65191895
                In [3]: dc.get_norm_eq('kb35')
                Out[3]: 1.92844639
                In [4]: dc.get_norm_eq('coss')
                Out[4]: 1.5
                In [5]: dc.get_norm_eq('rect')
                Out[5]: 1.0
                In [6]: dc.get_norm_eq('kbmd')
                Out[6]: 1.65994218
            Invalid inputs will prompt you for a new one.
                In [7]: dc.get_norm_eq('bob')
                Invalid window type. Please use one of the following:
                'rect', 'coss', 'kb25', 'kb35', or 'kbmd'
                Please enter a window type: kb25
                Out[7]: 1.65191895
        History:
            Translated from IDL: RJM - 2018/05/15 5:11 P.M.
    """
    if (type(wtype) == type('Hello')):
        wtype = wtype.replace(" ", "").lower()
        if (wtype in func_dict) and (np.size(wtype) == 1):
            norm_eq = func_dict[wtype]["normeq"]
        else:
            print("Invalid window type. Please use one of the following:")
            print("'rect', 'coss', 'kb25', 'kb35', or 'kbmd'")
            wtype   = input("Please enter a window type: ")
            wtype   = wtype.replace(" ", "").lower()
            wtype   = wtype.replace("'","")
            wtype   = wtype.replace('"',"")
            norm_eq = func_dict[wtype]["normeq"]
    else:
        print("Invalid window type. Please use one of the following:")
        print("'rect', 'coss', 'kb25', 'kb35', or 'kbmd'")
        wtype   = input("Please enter a window type: ")
        wtype   = wtype.replace(" ", "").lower()
        wtype   = wtype.replace("'","")
        wtype   = wtype.replace('"',"")
        norm_eq = func_dict[wtype]["normeq"]
    return norm_eq

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
    if not (type(T_in) in vtypes.atypes):
        sys.exit("Complex transmittance must be real or complex valued.")
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
    if not (type(T_in) in vtypes.atypes):
        sys.exit("Complex transmittance must be real or complex valued.")
    else:
        phase = np.arctan2(np.real(T_in),np.imag(T_in))
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
    if not (type(T_in) in vtypes.atypes):
        sys.exit("Complex transmittance must be real or complex.")
    elif (type(mu) not in vtypes.rtypes):
        sys.exit("mu must be real valued.")
    else:
        p           = power_func(T_in)
        crange      = (p>0).nonzero()
        tau         = np.zeros(np.size(p))
        tau[crange] = -mu[crange] * np.log(np.abs(p[crange]))
    return tau

def fresnel_cos(x_in):
    """
        Function:
            fresnel_cos
        Purpose:
            Compute the Fresnel cosine function.
        Variables:
            x_in:   A real or complex argument.
        Outputs:
            f_cos:  The fresnel cosine integral of x_in.
        Dependences:
            [1] diffcorr
            [2] numpy
            [3] sys
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
                is is useful for forward modeling problems in 
                radiative transfer and diffraction.
        References:
            [1] https://en.wikipedia.org/wiki/Fresnel_integral
            [2] https://en.wikipedia.org/wiki/Error_function
            [3] http://mathworld.wolfram.com/FresnelIntegrals.html
            [4] http://mathworld.wolfram.com/Erf.html
        Examples:
            Compute and plot the Fresnel Cosine integral.
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = dc.fresnel_cos(x)
                In [6]: plt.show(plt.plot(x,y))
        History:
            Translated from IDL:     RJM - 2018/05/14 2:13 P.M.
            Allow complex arguments: RJM - 2018/05/16 1:25 P.M.
    """
    ermes = "Input must be real or complex (Int, float, or array)"
    tx    = type(x_in)
    if not (tx in vtypes.atypes):
        sys.exit(ermes)
    f_cos = ((1-1j)/4.0)*erf((1+1j)*x_in*np.sqrt(np.pi)/2.0)+\
        ((1+1j)/4.0)*erf((1-1j)*x_in*np.sqrt(np.pi) / 2.0)
    if (not check_complex(x_in)):
        f_cos = np.real(f_cos)
    return f_cos

def fresnel_sin(x_in):
    """
        Function:
            fresnel_sin
        Purpose:
            Compute the Fresnel sine function.
        Variables:
            x_in:   A real or complex argument.
        Outputs:
            f_cos:  The fresnel sine integral of x_in.
        Dependences:
            [1] diffcorr
            [2] numpy
            [3] sys
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
            Compute and plot the Fresnel Cosine integral.
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: import matplotlib.pyplot as plt
                In [4]: x = np.array(range(0,10001))*0.01 - 50.0
                In [5]: y = dc.fresnel_sin(x)
                In [6]: plt.show(plt.plot(x,y))
        History:
            Translated from IDL: RJM - 2018/05/14 3:53 P.M.
            Allow complex arguments: RJM - 2018/05/16 1:26 P.M.
    """
    ermes = "freq_wav: Input must be a positive number (Int, float, or numpy)"
    tx    = type(x_in)
    if not (tx in vtypes.atypes):
        sys.exit(ermes)
    f_sin = ((1+1j)/4.0)*erf((1+1j)*x_in*np.sqrt(np.pi)/2.0)+\
    ((1-1j)/4.0)*erf((1-1j)*x_in*np.sqrt(np.pi)/2.0)
    if (not check_complex(x_in)):
        f_sin = np.real(f_sin)
    return f_sin

def sq_well_solve(x,a,b,F,Inverted=False):
    """
        Function:
            sq_well_solve
        Purpose:
            Computes the solution of diffraction through a square well.
        Variables:
            x: The independent variable.
            a: The LEFTMOST endpoint of the square well (Real number).
            b: The RIGHTMOST endpoint of the square well (Real number).
            F: The Fresnel scale.
        Output:
            H: Diffraction pattern of a square well on the interval [a,b].
        History:
            Translated from IDL: RJM - 2018/05/15 8:03 P.M.
    """
    if (not check_real(a)) or (not check_real(b)):
        sys.exit("Endpoints must be real valued")
    if (np.size(a) != 1) or (np.size(b) != 1):
        sys.exit("a and b must be a single values")
    if (not check_real(F)) or (np.size(F) != 1):
        sys.exit("F must be a single real value.")
    if (not check_real(x)):
        sys.exit("x must be real valued.")
    H = ((1 - 1j) / 2.) * (fresnel_cos((b - x) / F) - fresnel_cos((a - x) / F)\
    + 1j*(fresnel_sin((b - x) / F) - fresnel_sin((a - x) / F)))
    if not Inverted:
        H = 1-H
    return H

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
        wavfreq = spice.clight() / freqwav
    return wavfreq

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

def get_range_actual(rho,rng,w_vals):
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
    n_used      = finish - start
    return start, n_used

def fresnel_transform(T,ker,DX,F_SCALE):
    """
        Function: fresnel_transform
        Purpose:  Compute the approximate inverse of a Fresnel transform.
            Variables:
                T:   The input function (Diffraction data).
                KER: The Fresnel kernel.
                dx:  The size of a bin (r[1]-r[0]).
                F_SCALE: The Frensel scale. If not set, a default of 1 is used.
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
    if (not check_pos_real(F_SCALE)):
        sys.exit('F_SCALE must be positive')
    T_hat = np.sum(ker * T) * DX * (1.0-1.0j) / (2. * F_SCALE)
    return T_hat

def fresnel_inverse(T_hat,ker,dx,f_scale):
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
    T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
    return T

def psi_d1_phi_fast(r,r0,d,cb,cp,sp,cp0,sp0):
    """
        Function:
            psi_d1_phi_fast
        Purpose:
            Calculate dpsi/dphi from geometry variables using
            previously computed sines and cosines.
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
    if (not check_real(r)):
        sys.exit("r must be real valued")
    if (not check_real(r0)):
        sys.exit("r0 must be real valued")
    if (not check_real(d)):
        sys.exit("d must be real valued")
    if (not check_real(cb)):
        sys.exit("cos(b) must be real valued")
    if (not check_real(cp)):
        sys.exit("cos(phi) must be real valued")
    if (not check_real(sp)):
        sys.exit("sin(phi) must be real valued")
    if (not check_real(cp0)):
        sys.exit("cos(phi0) must be real valued")
    if (not check_real(sp0)):
        sys.exit("sin(phi0) must be real valued")
    xi   = (cb / d) * (r0*cp0 - r*cp)
    eta  = (r0**2 + r**2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / (d**2)
    v1   = r * cb * sp / d
    v2   = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d**2)
    psi_d1_phi_vals = (2.0*v1 + v2) / (2.0 * np.sqrt(1.0 + 2.0*xi + eta)) - v1
    return psi_d1_phi_vals

def psi_d1_phi(r,r0,d,b,phi,phi0):
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

def psi_d2_phi_fast(r,r0,d,cb,cp,sp,cp0,sp0):
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
            Translated from IDL: RJM - 2018/05/15 7:21 P.M.
    """
    if (not check_real(r)):
        sys.exit("r must be real valued")
    if (not check_real(r0)):
        sys.exit("r0 must be real valued")
    if (not check_real(d)):
        sys.exit("d must be real valued")
    if (not check_real(cb)):
        sys.exit("cos(b) must be real valued")
    if (not check_real(cp)):
        sys.exit("cos(phi) must be real valued")
    if (not check_real(sp)):
        sys.exit("sin(phi) must be real valued")
    if (not check_real(cp0)):
        sys.exit("cos(phi0) must be real valued")
    if (not check_real(sp0)):
        sys.exit("sin(phi0) must be real valued")
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

def psi_d2_phi(r,r0,d,b,phi,phi0):
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

def psi_factor_fast(r,r0,cb,cp0,sp0):
    """
        Function: psi_factor_fast
        Purpose:  Calculate the first iteration of Newton-Raphson for psi with
            respect to phi using previously calculated sines and cosines.
        Variables:
            r:   Ring radius variable.
            r0:  Ring intercept point.
            cb:  Cosine of ring opening angle.
            cp0: Cosine of ring azimuth angle.
            sp0: Sine of ring azimuth angle.
        History:
            Translated from IDL: RJM - Rough Draft - 2018/05/15 7:35 P.M.
    """
    if (not check_real(r)):
        sys.exit("r must be real valued")
    if (not check_real(r0)):
        sys.exit("r0 must be real valued")
    if (not check_real(cb)):
        sys.exit("cos(b) must be real valued")
    if (not check_real(cp0)):
        sys.exit("cos(phi0) must be real valued")
    if (not check_real(sp0)):
        sys.exit("sin(phi0) must be real valued")
    factor  = ((cb**2) * cp0 * sp0 / (1.0 - (cb**2) * (sp0**2))) * (r - r0) / r0
    return factor

def psi_factor(r,r0,b,phi0):
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
    if (not check_real(r)):
        sys.exit("r must be real valued")
    if (not check_real(r0)):
        sys.exit("r0 must be real valued")
    if (not check_real(b)):
        sys.exit("b must be real valued")
    if (not check_real(phi0)):
        sys.exit("phi0 must be real valued")
        cb      = np.cos(b)
        sp0     = np.sin(phi0)
        cp0     = np.cos(phi0)
        factor  = ((cb**2) * cp0 * sp0 / (1.0 - (cb**2) * (sp0**2))) * (r - r0) / r0
        return factor

def psi_fast(r,r0,d,cb,cp,sp,cp0,sp0):
    """
        Function: psi_fast
        Purpose:  Calculate psi from geometry variables.
        Variables:
            r:   Ring radius variable.
            r0:  Ring intercept point.
            D:   RIP-Spacecraft distance.
            cb:  Cosine of ring opening angle.
            cp:  Cosine of ring azimuth angle (variable).
            sp:  Sine of ring azimuth angle (variable).
            cp0: Cosine of ring azimuth angle (value).
            sp0: Sine of ring azimuth angle (value).
        History:
            Translated from IDL: RJM - 2018/05/15 7:48 P.M.
    """
    if (not check_real(r)):
        sys.exit("r must be real valued")
    if (not check_real(r0)):
        sys.exit("r0 must be real valued")
    if (not check_real(d)):
        sys.exit("d must be real valued")
    if (not check_real(cb)):
        sys.exit("cos(b) must be real valued")
    if (not check_real(cp)):
        sys.exit("cos(phi) must be real valued")
    if (not check_real(sp)):
        sys.exit("sin(phi) must be real valued")
    if (not check_real(cp0)):
        sys.exit("cos(phi0) must be real valued")
    if (not check_real(sp0)):
        sys.exit("sin(phi0) must be real valued")
    xi   = (cb / d) * (r0*cp0 - r*cp)
    eta  = ((r0**2) + (r**2) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d**2)
    psi_vals   = np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi)
    return psi_vals

def psi(r,r0,d,b,phi,phi0):
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
    cb   = np.cos(b)
    sp   = np.sin(phi)
    cp   = np.cos(phi)
    sp0  = np.sin(phi0)
    cp0  = np.cos(phi0)
    xi   = (cb / d) * (r0*cp0 - r*cp)
    eta  = ((r0**2) + (r**2) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d**2)
    psi_vals   = np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi)
    return psi_vals

def normalize(r,w_func,f_scale):
    """
        Function: normalize
        Purpose:  Compute the normalization factor used in the Fresnel
            Inversion process.
        Variables:
            r:       Rind radius.
            w_func:  The window (taper) function.
            f_scale: The Fresnel scale.
        History:
            Translated from IDL: RJM - 2018/05/15 8:21 P.M.
    """
    if (not check_real(r)):
        sys.exit("RHO must be real valued")
    if (not check_real(w_func)):
        sys.exit("W_FUNC must be real valued")
    if (not check_pos_real(f_scale)):
        sys.exit("F_SCALE must be a positive real number")
    if (np.size(r) != np.size(w_func)):
        sys.exit("RHO and W_FUNC have a different number of points")
    if (np.size(r) < 2.0):
        sys.exit("RHO needs to have at least 2 points")
    x         = r-np.mean(r)
    drho      = r[1]-r[0]
    psi       = (np.pi / 2.0) * ((x / f_scale)**2)
    ker       = np.exp(1j * psi)
    T1        = np.abs(np.sum(w_func * ker) * drho)
    norm_fact = np.sqrt(2.0) * f_scale / T1
    return norm_fact

def window_width(res,normeq,fsky,fres,rho_dot,sigma=False,bfac=True):
    """
        Function: window_width
        Purpose:  Compute the window width as a function of ring radius.
        Variables:
            res:     The requested resolution.
            normeq:  The normalized equivalent width. Unitless.
            fsky:    The sky frequency.
            fres:    The Fresnel scale.
            rho_dot: The time derivative of the ring radius.
        Output:
            w_vals:  The window width as a function of ring radius.
        History:
            Translated from IDL: RJM - 2018/05/15 8:38 P.M.
    """
    if bfac:
        if (not sigma):
            sigma = 2.e-13
        omega = 2.0 * np.pi * fsky
        alpha = (omega**2) * (sigma**2) / (2.0 * rho_dot)
        P     = res / (alpha * (fres**2))
        # The inverse exists only if P>1.
        if (np.min(P) < 1.0001):
            print("ERROR: Bad Points!")
            print((P < 1.0001).nonzero())
            print("Either rho_dot_kms_vals, F_km_vals, or res_km is to small.")
            print("Exclude points or request coarser resolution.")
            sys.exit("Illegal parameter values. Read error message.")
        # w_vals = np.abs(normeq * resolution_inverse(P) / alpha)
        w_vals = normeq*np.abs(resolution_inverse(P) / alpha)
    else:
        w_vals = 2.0*normeq*fres*fres/res
    return w_vals

def fresnel_forward_fast(rho_vals,F_vals,phi_rad_vals,b_rad_vals,d_vals,
    T_vals,lambda_vals,w_vals,dx,wtype,start,n_used,Normalize=True):
    """
        Procedure: fresnel_forward_fast
        Purpose:   Computes the forward model of diffraction from reconstructed
            data using a 'fast' method to speed up computation time. This is
            achieved by computing cosine and sine functions in the outer for
            loop, and then passing these computed values into the functions
            that need them.
        Variables:
            rho_vals:        Ring radius, in kilometers.
            F_vals:          Fresnel scale, in kilometers.
            phi_rad_vals:    Ring azimuth angle, in radians.
            b_rad_vals:      Ring opening angle, in radians.
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

    # Compute necessary variables.
    kD_vals   = 2. * np.pi * d_vals / lambda_vals
    w_max     = np.max(w_vals[start:start + n_used])
    nw_fwd    = np.ceil(w_max / (2. * dx))
    start_fwd = int(start + nw_fwd)
    n_fwd     = int(n_used - 2 * nw_fwd)
    cosb      = np.cos(b_rad_vals)
    cosphi0   = np.cos(phi_rad_vals)
    sinphi0   = np.sin(phi_rad_vals)

    # Compute forward model, point by point.
    T_hat_fwd_vals = T_vals * 0.0
    for i in range(n_fwd):
        center = start_fwd + i
        w      = w_vals[center]
        w_func = func_dict[wtype]["func"](w,dx)
        nw     = np.size(w_func)
        crange = np.array(range(int(center-(nw-1)/2),int(1+center+(nw-1)/2)))
        r      = rho_vals[crange]
        r0     = rho_vals[center]
        d      = d_vals[center]
        cb     = cosb[center]
        cp0    = cosphi0[center]
        sp0    = sinphi0[center]
        kD     = kD_vals[center]
        dphi_s_rad = psi_factor_fast(r,r0,cb,cp0,sp0)
        phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
        loop = 0

        # Perform Newton-Raphson on phi.
        while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
            cp         = np.cos(phi_s_rad)
            sp         = np.sin(phi_s_rad)
            psi_d1     = kD * psi_d1_phi_fast(r,r0,d,cb,cp,sp,cp0,sp0)
            psi_d2     = kD * psi_d2_phi_fast(r,r0,d,cb,cp,sp,cp0,sp0)
            dphi_s_rad = -psi_d1 / psi_d2
            phi_s_rad += dphi_s_rad
            loop      += 1
            if loop > 5:
                break

        # Compute psi and then compute the forward model.
        cp       = np.cos(phi_s_rad)
        sp       = np.sin(phi_s_rad)
        psi_vals = kD * psi_fast(r,r0,d,cb,cp,sp,cp0,sp0)
        ker      = wker(w_func,psi_vals)
        T        = T_vals[crange]
        F        = F_vals[center]
        T_hat_fwd_vals[center] = fresnel_transform(T,ker,dx,F)
        if Normalize:
            norm_factor = normalize(r,w_func,F)
            T_hat_fwd_vals[center] *= norm_factor
        print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Forward" \
        % (i,n_used,nw,loop),end="\r")
    print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Forward" \
        % (i,n_used,nw,loop))
    return T_hat_fwd_vals

def fresnel_forward(rho_vals,F_vals,phi_rad_vals,b_rad_vals,d_vals,
    T_vals,lambda_vals,w_vals,dx,wtype,start,n_used,Normalize=True):
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

    # Compute necessary variables.
    kD_vals   = 2. * np.pi * d_vals / lambda_vals
    w_max     = np.max(w_vals[start:start + n_used])
    nw_fwd    = np.ceil(w_max / (2. * dx))
    start_fwd = int(start + nw_fwd)
    n_fwd     = int(n_used - 2 * nw_fwd)

    # Compute forward model, point by point.
    T_hat_fwd_vals = T_vals * 0.0
    for i in range(n_fwd):
        center = start_fwd + i
        w      = w_vals[center]
        w_func = func_dict[wtype]["func"](w,dx)
        nw     = np.size(w_func)
        crange = np.array(range(int(center-(nw-1)/2),int(1+center+(nw-1)/2)))
        r      = rho_vals[crange]
        r0     = rho_vals[center]
        d      = d_vals[center]
        b      = b_rad_vals[center]
        phi0   = phi_rad_vals[center]
        kD     = kD_vals[center]
        dphi_s_rad = psi_factor(r,r0,b,phi0)
        phi_s_rad  = phi0 - dphi_s_rad
        loop = 0

        # Perform Newton-Raphson on phi.
        while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
            psi_d1     = kD * psi_d1_phi(r,r0,d,b,phi_s_rad,phi0)
            psi_d2     = kD * psi_d2_phi(r,r0,d,b,phi_s_rad,phi0)
            dphi_s_rad = -psi_d1 / psi_d2
            phi_s_rad += dphi_s_rad
            loop      += 1
            if loop > 5:
                break

        # Compute psi and then compute the forward model.
        psi_vals = kD * psi(r,r0,d,b,phi_s_rad,phi0)
        ker      = wker(w_func,psi_vals)
        T        = T_vals[crange]
        F        = F_vals[center]
        T_hat_fwd_vals[center] = fresnel_transform(T,ker,dx,F)
        if Normalize:
            norm_factor = normalize(r,w_func,F)
            T_hat_fwd_vals[center] *= norm_factor
        print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Normal Forward" \
        % (i,n_used,nw,loop),end="\r")
    print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Normal Forward" \
    % (i,n_used,nw,loop))
    return T_hat_fwd_vals

def fresnel_inversion_fast(rho_vals,F_vals,phi_rad_vals,b_rad_vals,d_vals,
    T_hat_vals,lambda_vals,w_vals,dx,wtype,start,n_used,Normalize=True):
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
    kD_vals   = 2. * np.pi * d_vals / lambda_vals
    cosb      = np.cos(b_rad_vals)
    cosphi0   = np.cos(phi_rad_vals)
    sinphi0   = np.sin(phi_rad_vals)
    # Define functions
    fw        = func_dict[wtype]["func"]
    psid1     = psi_d1_phi_fast
    psid2     = psi_d2_phi_fast
    psifac    = psi_factor_fast
    finv      = fresnel_inverse
    psif      = psi_fast
    nrm       = normalize
    # Calculate the corrected complex amplitude, point by point
    T_vals     = T_hat_vals * 0.0
    nw1        = 0
    phi_s_rad1 = None
    for i in np.arange(n_used+1):
        center = start+i
        w      = w_vals[center]
        w_func = fw(w,dx)
        nw     = np.size(w_func)
        crange = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
        r      = rho_vals[crange]
        r0     = rho_vals[center]
        d      = d_vals[center]
        cb     = cosb[center]
        cp0    = cosphi0[center]
        sp0    = sinphi0[center]
        kD     = kD_vals[center]

        if (nw == nw1):
            phi_s_rad  = phi_s_rad1
            cp         = np.cos(phi_s_rad)
            sp         = np.sin(phi_s_rad)
            psi_d1     = psid1(r,r0,d,cb,cp,sp,cp0,sp0)
            psi_d2     = psid2(r,r0,d,cb,cp,sp,cp0,sp0)
            dphi_s_rad = -psi_d1 / psi_d2
        else:
            dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
            phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
        phi_s_rad1 = phi_s_rad
        nw1        = nw
        loop = 0

        # Perform Newton-Raphson on phi.
        while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
            cp         = np.cos(phi_s_rad)
            sp         = np.sin(phi_s_rad)
            psi_d1     = psid1(r,r0,d,cb,cp,sp,cp0,sp0)
            psi_d2     = psid2(r,r0,d,cb,cp,sp,cp0,sp0)
            dphi_s_rad = -psi_d1 / psi_d2
            phi_s_rad += dphi_s_rad
            loop      += 1
            if loop > 5:
                break
        
        # Compute psi and then compute the forward model.
        cp       = np.cos(phi_s_rad)
        sp       = np.sin(phi_s_rad)
        psi_vals = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
        ker      = wker(w_func,-psi_vals)
        T_hat    = T_hat_vals[crange]
        F        = F_vals[center]
        T_vals[center] = finv(T_hat,ker,dx,F)
        if Normalize:
            T_vals[center] *= nrm(r,w_func,F)
        #print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
        #% (i,n_used,nw,loop),end="\r")
    #print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
    #% (i,n_used,nw,loop))
    return T_vals

def fresnel_inversion(rho_vals, F_vals, phi_rad_vals, b_rad_vals,d_vals,
    T_hat_vals, lambda_vals, w_vals, dx, wtype, start, n_used, Normalize=True):
    """
        Function:
            fresnel_inversion 
        Purpose:
            Computes the fresnel inversion from a set of 
            diffracted data using a 'fast' method to speed up 
            computation time. This is achieved by computing cosine
            and sine function in the outer for loop, and then passing
            these computed values into the functions that need them.
            The normal version passes the arguments to the functions,
            and then cosines and sines are computed within the
            function. For small data sets or coarser resolutions,
            the normal version is faster. Normal and fast have
            identical outputs.
        Variables:
            rho_vals:           Ring radius, in kilometers.
            F_vals:             Fresnel scale, in kilometers.
            phi_rad_vals:       Ring azimuth angle, in radians.
            B_rad_vals:         Ring opening angle, in radians.
            lambda_sky_vals:    Wavelength of recieved signal, in
                                kilometers.
            D_vals:             Spacecraft-RIP distance, in
                                kilometers.
            dx:                 Sampling spacing, in kilometers.
            T_vals:             Reconstructed complex transmittance.
            w_vals:             Window width, in kilometers.
            wtype:              Window used in reconstruction,
                                string.
            start:              Starting point of reconstructed data.
            n_used:             Number of reconstructed points.
        Keywords
            Normalize:  Parameter for normalizing the complex
                        transmittance by the window function that
                        is used. Default is True. Set to False to
                        skip this feature.
        Output:
            T_vals:     Reconstructed Complex transmittance.
        History:
            Translated from IDL: RJM - 2018/05/16 6:52 A.M.
    """
    # Compute necessary variables.
    kD_vals   = 2. * np.pi * d_vals / lambda_vals
 
    # Calculate the corrected complex amplitude, point by point
    T_vals = T_hat_vals * 0.0
    for i in range(n_used):
        center = start + i
        w      = w_vals[center]
        w_func = func_dict[wtype]["func"](w,dx)
        nw     = np.size(w_func)
        crange = np.array(range(int(center-(nw-1)/2),int(1+center+(nw-1)/2)))
        r      = rho_vals[crange]
        r0     = rho_vals[center]
        d      = d_vals[center]
        b      = b_rad_vals[center]
        phi0   = phi_rad_vals[center]
        kD     = kD_vals[center]
        dphi_s_rad = psi_factor(r,r0,b,phi0)
        phi_s_rad  = phi0 - dphi_s_rad
        loop = 0

        # Perform Newton-Raphson on phi.
        while (np.max(np.abs(dphi_s_rad)) > 1.e-10):
            psi_d1     = kD * psi_d1_phi(r,r0,d,b,phi_s_rad,phi0)
            psi_d2     = kD * psi_d2_phi(r,r0,d,b,phi_s_rad,phi0)
            dphi_s_rad = -psi_d1 / psi_d2
            phi_s_rad += dphi_s_rad
            loop      += 1
            if loop > 5:
                break
        
        # Compute psi and then compute the forward model.
        psi_vals = kD * psi(r,r0,d,b,phi_s_rad,phi0)
        ker      = wker(w_func,-psi_vals)
        T        = T_hat_vals[crange]
        F        = F_vals[center]
        T_vals[center] = fresnel_inverse(T,ker,dx,F)
        if Normalize:
            norm_factor = normalize(r,w_func,F)
            T_vals[center] *= norm_factor
        print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Normal Inversion" \
        % (i,n_used,nw,loop),end="\r")
    print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Normal Inversion" \
    % (i,n_used,nw,loop))
    return T_vals