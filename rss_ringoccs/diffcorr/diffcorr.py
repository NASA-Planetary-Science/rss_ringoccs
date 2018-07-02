"""
    Module Name:
        diffcorr
    Purpose:
        Provide functions and classes that aid in the process of
        Diffraction Correction / Fresnel Inversion. Additional
        functions for the purpose of forward modelling of
        reconstructed data and diffraction modelling are included.
        Several low-level functions that perform error checks for
        the main functions also exist, as well as functions that
        enable running shell scripts in Python.
    
    Window (Taper) Functions:
        rect................Rectangular window.
        coss................Squared cossine window.
        kb20................Kaiser-Bessel 2.0 window.
        kb25................Kaiser-Bessel 2.5 window.
        kb35................Kaiser-Bessel 3.5 window.
        kbmd20..............Modified Kaiser-Bessel 2.0 window.
        kbmd25..............Modified Kaiser-Bessel 2.5 window.

    Error Checks:
        check_boole.........Checks for Boolean input.
        check_ternary.......Checks for Ternary input.
        check_pos_real......Checks for a single positive real input.
        check_real..........Checks if input is real (Number/Array).
        check_complex.......Checks if input is complex.
    
    Special Functions:
        fresnel_sin.........The Fresnel sine integral.
        fresnel_cos.........The Fresnel cosine integral.
        sq_well_solve.......Diffraction pattern through square well.

    Mathematical Functions:
        compute_norm_eq.....Computes the normalized equivalent width.
        get_norm_eq.........Quickly retrieve pre-computed normalized
                            equivalent widths from strings with the
                            name of common window functions.
        resolution_inverse..Computes the inverse of the function
                            y = x/(exp(-x)+x-1)
        power_func..........Compute power from complex transmittance.
        phase_func..........Compute phase from complex transmittance.
        tau_func............Compute normalized optical depth from the
                            complex transmittance.
        wker................Computes a weighted kernel function.
        freq_wav............Convert frequency to wavelength, and
                            vice-versa. Kilometers or Herts only.
        fresnel_scale.......Compute the Fresnel scale.
    
    Miscellaneous Functions:
        get_range_request...Computes an array of the form [a,b] from
                            a given array, list, or from a set of
                            allowed strings.
        get_range_actual....Given an array of numbers (usually the
                            radial range of the data), a range
                            request, and a window width, compute the
                            allowed range of processing.
"""

# Import dependencies for the diffcorrpy module
import time
import os
import sys
import platform
import numpy as np
from scipy.special import lambertw, iv, erf
from scipy import interpolate
from scipy.constants import speed_of_light
import pdb

sys.path.append("../../")
from rss_ringoccs.tools import check_boole,check_ternary
from rss_ringoccs.tools import check_real,check_pos_real
from rss_ringoccs.tools import check_complex,get_geo

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
    if not check_real(w_func):
        raise TypeError("Input must be real valued")
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
    if (not check_real(x)) and (not check_complex(x)):
        raise TypeError("Input must be real or complex")
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
    if (not check_real(rngreq)) and (ti != type([1])) and (ti != type("Hi")):
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
            print("'rect', 'coss', 'kb25', 'kb20', kb35', or 'kbmd'")
            wtype   = input("Please enter a window type: ")
            wtype   = wtype.replace(" ", "").lower()
            wtype   = wtype.replace("'","")
            wtype   = wtype.replace('"',"")
            norm_eq = func_dict[wtype]["normeq"]
    else:
        print("Invalid window type. Please use one of the following:")
        print("'rect', 'coss', 'kb25', 'kb20', kb35', or 'kbmd'")
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
    if (not check_real(x_in)) and (not check_complex(x_in)):
        raise TypeError("Input must be real or complex.")
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
    if (not check_real(x_in)) and (not check_complex(x_in)):
        raise TypeError("Input must be real or complex.")
    f_sin = ((1+1j)/4.0)*erf((1+1j)*x_in*np.sqrt(np.pi)/2.0)+\
    ((1-1j)/4.0)*erf((1-1j)*x_in*np.sqrt(np.pi)/2.0)
    if (not check_complex(x_in)):
        f_sin = np.real(f_sin)
    return f_sin

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
        raise ValueError("Input must be two positive real numbers")
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
    x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
    w_func = np.cos(np.pi * x / w_in)**2
    return w_func

def kb20(w_in, dx):
    nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
    x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
    alpha  = 2.0*np.pi
    w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
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
    x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
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
    x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
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
    x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
    alpha  =  al * np.pi
    w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
    return w_func

def kbmd20(w_in, dx):
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
    alpha  = 2.0 * np.pi
    w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
    return w_func

def kbmd25(w_in, dx):
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
    n_used      = 1 + (finish - start)
    return start, n_used

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

def psi_d1_phi(object):
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

def fresnel_inverse_fft(T_hat,ker,dx,f_scale):
    fft_t_hat       = np.fft.fft(T_hat)
    fft_conv        = np.fft.fft(ker)
    inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
    inv_t_hat      *= dx*(np.complex(1.0,1.0))/(2.0*f_scale)
    T               = inv_t_hat[int((nw-1)/2)+1]
    return T

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
    eta  = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d*d)
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

func_dict = {
    "rect"      : {"func" : rect, "normeq" : 1.00000000},
    "coss"      : {"func" : coss, "normeq" : 1.50000000},
    "kb20"      : {"func" : kb20, "normeq" : 1.49634231},
    "kb25"      : {"func" : kb25, "normeq" : 1.65191895},
    "kb35"      : {"func" : kb35, "normeq" : 1.92844639},
    "kbmd20"    : {"func" : kbmd20, "normeq" : 1.52048174},
    "kbmd25"    : {"func" : kbmd25, "normeq" : 1.65994218}
    }

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
        w_vals = normeq*np.abs(resolution_inverse(P) / alpha)
    else:
        w_vals = 2.0*normeq*fres*fres/res
    return w_vals

def fresnel_forward_fast(rho_vals,F_vals,phi_rad_vals,B_rad_vals,d_vals,
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
    cosb      = np.cos(B_rad_vals)
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

def fresnel_forward(rho_vals,F_vals,phi_rad_vals,B_rad_vals,d_vals,
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
        b      = B_rad_vals[center]
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

def fresnel_inversion_fast(rho_vals,F_vals,phi_rad_vals,B_rad_vals,d_vals,
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
    cosb      = np.cos(B_rad_vals)
    cosphi0   = np.cos(phi_rad_vals)
    sinphi0   = np.sin(phi_rad_vals)
    # Define functions
    fw        = __func_dict[wtype]["func"]
    psid1     = __psid1fast
    psid2     = __psid2fast
    psifac    = __psifacfast
    finv      = __fresinv
    psif      = __psifast
    nrm       = __normalize
    # Calculate the corrected complex amplitude, point by point
    T_vals = T_hat_vals * 0.0
    nw1    = 0
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

def fresnel_inversion(rho_vals, F_vals, phi_rad_vals, B_rad_vals,d_vals,
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
        b      = B_rad_vals[center]
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

class rec_data(object):
    """
        Class: rec_data
        Purpose: This object contains all of input and output variables from
        diffraction reconstruction, including geometry, data, and reconstructed
        data.
        Attributes:
            rho_km_vals (np.ndarray): Ring intercept radius values, in km,
                at final spacing specified in dr_km
            p_norm_vals (np.ndarray): Power normalized to 1. This is the diffraction
                pattern that is input to the Fresnel inversion step
            phase_rad_vals (np.ndarray): Phase of the complex signal, in radians.
                This is the other part of the diffraction pattern that is input to
                the Fresnel Inversion step
            B_rad_vals (np.ndarray): Ring opening angle of Saturn
            D_km_vals (np.ndarray): Spacecraft-RIP (Ring Intercept Point) distance
            f_sky_hz_vals (np.ndarray): Predicted sky frequency, in Hz.
            phi_rad_vals (np.ndarray): Observed ring azimuth, in radians
            rho_dot_kms_vals (np.ndarray): Ring intercept point velocity
	"""

    def __init__(self, NormDiff,res,wtype,bfac=True):
        """
            Arguments:
                NormDiff: Instance of the class NormDiff containing 
        """

        self.res                     = None
        self.wtype                   = None
        self.rho_km_vals             = None
        self.p_norm_vals             = None
        self.phase_rad_vals          = None
        self.B_rad_vals              = None
        self.D_km_vals               = None
        self.f_sky_hz_vals           = None
        self.phi_rad_vals            = None
        self.rho_dot_kms_vals        = None
        self.T_hat_vals              = None
        self.F_km_vals               = None
        self.w_km_vals               = None
        self.mu_vals                 = None
        self.lambda_sky_km_vals      = None
        self.dx_km                   = None
        self.norm_eq                 = None
        self.t_oet_spm_vals          = None
        self.t_ret_spm_vals          = None
        self.t_set_spm_vals          = None
        self.rho_corr_pole_km_vals   = None
        self.rho_corr_timing_km_vals = None
        self.phi_rl_rad_vals         = None

        self.res   = res
        self.wtype = wtype.replace(" ", "").lower()

        self.__set_attributes(NormDiff)
        self.__get_occ_type()

        n_rho = np.size(self.rho_km_vals)
        self.__error_check(n_rho)

        self.lambda_sky_km_vals = freq_wav(self.f_sky_hz_vals)
        self.dx_km              = self.rho_km_vals[1] - self.rho_km_vals[0]
        self.T_hat_vals         = np.sqrt(np.abs(self.p_norm_vals)) * \
            np.exp(1j * self.phase_rad_vals)
        self.F_km_vals          = fresnel_scale(self.lambda_sky_km_vals,
            self.D_km_vals,self.phi_rad_vals,self.B_rad_vals)
        self.norm_eq            = get_norm_eq(self.wtype)
        self.w_km_vals          = window_width(self.res,self.norm_eq,
            self.f_sky_hz_vals,self.F_km_vals,self.rho_dot_kms_vals,bfac=bfac)
        self.mu_vals            = np.sin(np.abs(self.B_rad_vals))
    
    def __set_attributes(self, NormDiff):
        self.rho_km_vals      = NormDiff.rho_km_vals
        self.p_norm_vals      = NormDiff.p_norm_vals
        self.phase_rad_vals   = -NormDiff.phase_rad_vals
        self.B_rad_vals		  = NormDiff.B_rad_vals
        self.D_km_vals		  = NormDiff.D_km_vals
        self.f_sky_hz_vals    = NormDiff.f_sky_hz_vals
        self.phi_rad_vals     = NormDiff.phi_rad_vals
        self.rho_dot_kms_vals = NormDiff.rho_dot_kms_vals
        try:
            self.t_oet_spm_vals          = NormDiff.t_oet_spm_vals         
            self.t_ret_spm_vals          = NormDiff.t_ret_spm_vals         
            self.t_set_spm_vals          = NormDiff.t_set_spm_vals         
            self.rho_corr_pole_km_vals   = NormDiff.rho_corr_pole_km_vals
            self.rho_corr_timing_km_vals = NormDiff.rho_corr_timing_km_vals
            self.phi_rl_rad_vals         = NormDiff.phi_rl_rad_vals
        except AttributeError:
            try:
                self.t_oet_spm_vals          = NormDiff.spm_vals
                self.t_ret_spm_vals          = NormDiff.t_ret_spm_vals         
                self.t_set_spm_vals          = NormDiff.t_set_spm_vals         
                self.rho_corr_pole_km_vals   = NormDiff.rho_corr_pole_km_vals
                self.rho_corr_timing_km_vals = NormDiff.rho_corr_timing_km_vals
                self.phi_rl_rad_vals         = NormDiff.phi_rl_rad_vals
            except AttributeError:
                pass

    
    def __error_check(self,N_RHO):
        error_code = 0
        if (np.size(self.phase_rad_vals) != N_RHO):
            print("Bad Input: len(phase_rad_vals) != len(rho_km_vals)")
            error_code += 1
        if (np.size(self.p_norm_vals) != N_RHO):
            print("Bad Input: len(p_norm_vals) != len(rho_km_vals)")
            error_code += 2
        if (np.size(self.phi_rad_vals) != N_RHO):
            print("Bad Input: len(phi_rad_vals) != len(rho_km_vals)")
            error_code += 4
        if (np.size(self.B_rad_vals) != N_RHO):
            print("Bad Input: len(B_rad_vals) != (rho_km_vals)")
            error_code += 8
        if (np.size(self.f_sky_hz_vals) != N_RHO):
            print("Bad Input: len(f_sky_Hz_vals) != len(rho_km_vals)")
            error_code += 16
        if (np.size(self.D_km_vals) != N_RHO):
            print("Bad Input: len(D_km_vals) != len(rho_km_vals)")
            error_code += 32
        if (np.size(self.rho_dot_kms_vals) != N_RHO):
            print("Bad Input: len(D_km_vals) != len(rho_km_vals)")
            error_code += 64
        if (error_code != 0):
            raise ValueError("Bad input. Read error message.")

    def __get_occ_type(self):
        drho = [np.min(self.rho_dot_kms_vals),np.max(self.rho_dot_kms_vals)]
        dx   = self.rho_km_vals[1]-self.rho_km_vals[0]
        if (drho[0] < 0) and (drho[1] > 0):
            raise ValueError("drho/dt had positive and negative values.")
        if (dx > 0) and (drho[1] < 0):
            self.rho_dot_kms_vals=np.abs(self.rho_dot_kms_vals)
        if (dx < 0):
            self.rho_km_vals      = self.rho_km_vals[::-1]
            self.phase_rad_vals   = self.phase_rad_vals[::-1]
            self.p_norm_vals      = self.p_norm_vals[::-1]
            self.phi_rad_vals     = self.phi_rad_vals[::-1]
            self.B_rad_vals       = self.B_rad_vals[::-1]
            self.f_sky_hz_vals    = self.f_sky_hz_vals[::-1]
            self.D_km_vals        = self.D_km_vals[::-1]
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])

class diffraction_correction(object):
    """
        Class:
            diffraction_correction
        Purpose:
            Perform diffraction correction for a ring occultation
            on a data set that is a near radially symmetric function
            of the ring radius, or ring intercept point (RIP).
        Arguments:
            dat:    
                The data set, usually an instance of the NormDiff
                class from the rss_ringoccs Calibration subpackage.
                This instance MUST contain the following attributes
                and MUST have the same names.
                    rho_km_vals:      Ring Radius (km)
                    phi_rad_vals:     Ring Azimuth Angle (Radians)
                    p_norm_vals:      Normalized Power
                    phase_rad_vals:   Phase (Radians)
                    B_rad_vals:       Elevation Angle (Radians)
                    D_km_vals:        RIP-Distance (km)
                    f_sky_hz_vals:    Sky Frequency (Hertz)
                    rho_dot_kms_vals: RIP-velocity (km/s)
                    history:          History dictionary
            res:    
                The requested resolution for processing (km). This
                must be a positive real number, that is, a positive
                floating point number or integer.
        Keywords:
            rng:    
                The request range for diffraction correction.
                Preferred input is rng = [a,b]. Arrays are
                allowed and the range will be set as:
                    rng = [MIN(array),MAX(array)]
                Finally, certain strings containing a few of the
                regions of interests within the rings of Saturn
                are allowed. Permissible strings are:
                    'maxwell', 'titan', 'huygens', and 'encke'.
                Strings are neither case nor space sensitive.
                For other planets use rng = [a,b]. Default value
                is set to 'all' which processing [65,000,140,000]
                Values MUST be set in kilometers.
            wtype:  
                The requested tapering function for diffraction
                correction. A string with several allowed inputs:
                    'rect'      Rectangular Window.
                    'coss'      Squares Cosine Window.
                    'kb20'      Kaiser-Bessel 2.0 Window.
                    'kb25'      Kaiser-Bessel 2.5 Window.
                    'kb35'      Kaiser-Bessel 3.5 Window.
                    'kbmd20'    Modified kb20 Window.
                    'kbmd25'    Modified kb25 Window.
                The variable is neither case nor space sensitive.
                Default window is set to 'kb25'
            fwd:    
                A Boolean for determining whether or not
                forward modelling will be computed. This is good
                starting point for deciding if the diffraction
                correction is physically significant or valid. If
                the reconstruction is good, the forward model
                should reproduce the p_norm_vals attribute from
                the input dat instance. Default is set to False.
            norm:
                A Boolean for determining whether or not the
                reconstructed complex transmittance is normalize
                by the window width. This normalization is the
                complex transmittance that is computed by using
                free space divided by the complex transmittance
                that is computed using free space weighted by the
                selected tapering function. Default is True.
            bfac:
                A Boolean for determining whether or not the
                'b' factor in the window width computation is
                used. This is equivalent to setting the Allen
                Deviation from the spacecraft to a positive value
                or to zero. If set to False, the Allen Deviation
                is assumed to be zero. If set to True the Allen
                Deviation is set to 2e-13. If set to a positive
                real number, the Allen Deviation will be assumed
                to be that real number. Default is True.
            fft:    
                A Boolean for determining whether or not FFT's will
                be used for computing the complex transmittance. The
                use of FFT's assumes that the geometry of the system
                is such that the integral that is used to compute the
                complex transmittance is of the form of a
                convolution, and that the convolution theorem may be
                applied to it. Default is set to False.
            psitype:
                A string for determining what approximation to the
                geometrical 'psi' function is used. Several strings
                are allowed:
                    'full'      No Approximation is applied.
                    'taylor2'   Second order Taylor Series.
                    'taylor3'   Third order Taylor Series.
                    'taylor4'   Fourth order Taylor Series.
                    'MTR2'      Second Order Series from MTR86.
                    'MTR3'      Third Order Series from MTR86.
                    'MTR4'      Fourth Order Series from MTR86.
                The variable is neither case nor space sensitive.
                Default is set to 'full'.
            verbose:
                A Boolean for determining if various pieces of
                information are printed to the screen or not.
                Default is False.
        Outputs:
            T_hat_vals:
                Complex transmittance of the diffracted data.
            F_km_vals:
                Fresnel scale (km).
            w_km_vals:
                Window width as a function of ring radius (km).
            mu_vals:
                The sine of the elevation angle.
            lambda_sky_km_vals:
                Wavelength of the recieved signal (km).
            dx_km:  
                Radial spacing between points (km).
            norm_eq:        
                Normalized equivalent width of the window function.
            n_used:
                Number of points that were processed (integer).
            start:
                Starting point used for processing (integer).
            T_vals:
                Complex transmittance of reconstructed data.
            power_vals:
                Normalized power of the reconstructed data.
            tau_vals:
                Normalized optical depth of the reconstructed data.
            phase_vals:
                Phase of the reconstructed data (Radians).
            p_norm_fwd_vals:
                Normalized power of forward model (fwd=True).
            T_hat_fwd_vals:
                Complex transmittance of forward model (fwd=True).
            phase_fwd_vals:
                Phase of forward model (fwd=True).
            history:
                History dictionary of the runtime OS settings.
        Dependencies:
            [1] numpy
            [2] scipy
            [3] diffcorr
            [4] 
        Notes:
            [1] 
        References:

        Examples:

        History:
            Created: RJM - 2018/05/16 5:40 P.M.
    """
    def __init__(self,dat,res,rng="all",wtype="kb25",fwd=False,
        norm=True,bfac=True,fft=False,psitype="full",verbose=True):
        t1       = time.time()

        self.res                        = None
        self.wtype                      = None
        self.rng                        = None
        self.rho_km_vals                = None
        self.p_norm_vals                = None
        self.phase_rad_vals             = None
        self.B_rad_vals                 = None
        self.D_km_vals                  = None
        self.f_sky_hz_vals              = None
        self.phi_rad_vals               = None
        self.rho_dot_kms_vals           = None
        self.T_hat_vals                 = None
        self.F_km_vals                  = None
        self.w_km_vals                  = None
        self.mu_vals                    = None
        self.lambda_sky_km_vals         = None
        self.dx_km                      = None
        self.norm_eq                    = None
        self.n_used                     = None
        self.start                      = None
        self.T_vals                     = None
        self.power_vals                 = None
        self.tau_vals                   = None
        self.phase_vals                 = None
        self.p_norm_fwd_vals            = None
        self.T_hat_fwd_vals             = None
        self.phase_fwd_vals             = None
        self.norm                       = None
        self.fwd                        = None
        self.fft                        = None
        self.bfac                       = None
        self.psitype                    = None
        self.history                    = None
        self.dathist                    = None
        self.tau_threshold_vals         = None
        self.t_oet_spm_vals             = None
        self.t_ret_spm_vals             = None
        self.t_set_spm_vals             = None
        self.rho_corr_pole_km_vals      = None
        self.rho_corr_timing_km_vals    = None
        self.phi_rl_rad_vals            = None

        if not check_boole(fwd):
            raise TypeError("fwd must be Boolean: True/False")
        if not check_boole(norm):
            raise TypeError("norm must be Boolean: True/False")
        if not check_boole(bfac):
            raise TypeError("bfac must be Boolean: True/False")
        if not check_boole(fft):
            raise TypeError("fft must be Boolean: True/False")
        if not check_boole(verbose):
            raise TypeError("verbose must be Boolean: True/False")
        if (type(wtype) != type("Hi")):
            raise TypeError("wtype must be a string: 'coss'")

        recdata                 = rec_data(dat,res,wtype,bfac=bfac)
        self.res                = res
        self.wtype              = wtype
        self.rng                = rng
        self.norm               = norm
        self.fwd                = fwd
        self.fft                = fft
        self.bfac               = bfac
        self.psitype            = psitype
        self.dathist            = dat.history
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.B_rad_vals         = recdata.B_rad_vals
        self.D_km_vals          = recdata.D_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq

        self.t_oet_spm_vals          = recdata.t_oet_spm_vals         
        self.t_ret_spm_vals          = recdata.t_ret_spm_vals         
        self.t_set_spm_vals          = recdata.t_set_spm_vals         
        self.rho_corr_pole_km_vals   = recdata.rho_corr_pole_km_vals     
        self.rho_corr_timing_km_vals = recdata.rho_corr_timing_km_vals
        self.phi_rl_rad_vals         = recdata.phi_rl_rad_vals

        del recdata,res,wtype,rng,norm,fwd,fft,bfac,psitype

        self.rng                = get_range_request(self.rng)
        self.start,self.n_used  = get_range_actual(self.rho_km_vals,
            self.rng,self.w_km_vals)

        self.__trim_inputs()
        if verbose:
            if (self.psitype == "full"):
                if self.norm: self.T_vals = self.__finv_n_v()
                else:    self.T_vals = self.__finv_v()
                if self.fwd:
                    if self.norm: self.T_hat_fwd_vals = self.__ffwd_n_v()
                    else:    self.T_hat_fwd_vals = self.__ffwd_v()
                    self.p_norm_fwd_vals = power_func(self.T_hat_fwd_vals)
                    self.phase_fwd_vals  = phase_func(self.T_hat_fwd_vals)
            else:
                if self.norm: self.T_vals = self.__finv_n_p_v(self.psitype)
                else:    self.T_vals = self.__finv_p_v(self.psitype)
                if self.fwd:
                    if self.norm: self.T_hat_fwd_vals = self.__ffwd_n_v()
                    else:    self.T_hat_fwd_vals = self.__ffwd_v()
                    self.p_norm_fwd_vals = power_func(self.T_hat_fwd_vals)
                    self.phase_fwd_vals  = phase_func(self.T_hat_fwd_vals)
        else:
            if (self.psitype == "full"):
                if self.norm: self.T_vals = self.__finv_n()
                else:    self.T_vals = self.__finv()
                if self.fwd:
                    if self.norm: self.T_hat_fwd_vals = self.__ffwd_n()
                    else:    self.T_hat_fwd_vals = self.__ffwd()
                    self.p_norm_fwd_vals = power_func(self.T_hat_fwd_vals)
                    self.phase_fwd_vals  = phase_func(self.T_hat_fwd_vals)
            else:
                if self.norm: self.T_vals = self.__finv_n_p(self.psitype)
                else:    self.T_vals = self.__finv_p(self.psitype)
                if self.fwd:
                    if self.norm: self.T_hat_fwd_vals = self.__ffwd_n()
                    else:    self.T_hat_fwd_vals = self.__ffwd()
                    self.p_norm_fwd_vals = power_func(self.T_hat_fwd_vals)
                    self.phase_fwd_vals  = phase_func(self.T_hat_fwd_vals)

        self.power_vals = power_func(self.T_vals)
        self.phase_vals = -phase_func(self.T_vals)
        self.tau_vals   = tau_func(self.T_vals,self.mu_vals)
        self.__trim_attributes(self.fwd)
        del verbose

        self.tau_threshold_vals = np.zeros(np.size(self.rho_km_vals))

        self.history = self.__write_tau_hist()
        t2 = time.time()
        sys.stdout.write("\033[K")
        print("Computation Time: ",t2-t1,end="\r")
    
    def __trim_attributes(self,fwd):
        start  = self.start
        n_used = self.n_used
        crange = np.arange(n_used)+start

        self.rho_km_vals             = self.rho_km_vals[crange]
        self.p_norm_vals             = self.p_norm_vals[crange]
        self.phase_rad_vals          = self.phase_rad_vals[crange]
        self.B_rad_vals              = self.B_rad_vals[crange]
        self.D_km_vals               = self.D_km_vals[crange]
        self.f_sky_hz_vals           = self.f_sky_hz_vals[crange]
        self.phi_rad_vals            = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals        = self.rho_dot_kms_vals[crange]
        self.T_hat_vals              = self.T_hat_vals[crange]
        self.F_km_vals               = self.F_km_vals[crange]
        self.w_km_vals               = self.w_km_vals[crange]
        self.mu_vals                 = self.mu_vals[crange]
        self.lambda_sky_km_vals      = self.lambda_sky_km_vals[crange]
        self.T_vals                  = self.T_vals[crange]
        self.power_vals              = self.power_vals[crange]
        self.tau_vals                = self.tau_vals[crange]
        self.phase_vals              = self.phase_vals[crange]
        try:
            self.t_oet_spm_vals          = self.t_oet_spm_vals[crange]
            self.t_ret_spm_vals          = self.t_ret_spm_vals[crange]
            self.t_set_spm_vals          = self.t_set_spm_vals[crange] 
            self.rho_corr_pole_km_vals   = self.rho_corr_pole_km_vals[crange]    
            self.rho_corr_timing_km_vals = self.rho_corr_timing_km_vals[crange]
            self.phi_rl_rad_vals         = self.phi_rl_rad_vals[crange]
        except TypeError:
            pass
        if fwd:
            self.p_norm_fwd_vals = self.p_norm_fwd_vals[crange]
            self.T_hat_fwd_vals  = self.T_hat_fwd_vals[crange]
            self.phase_fwd_vals  = self.phase_fwd_vals[crange]
    
    def __trim_inputs(self):
        start   = self.start
        n_used  = self.n_used
        rho     = self.rho_km_vals
        rstart  = rho[start]
        rfin    = rho[start+n_used+1]
        w_vals  = self.w_km_vals
        w       = np.ceil(np.max(w_vals[start:start+n_used+1])/2.0)
        nst     = np.min((rho>=(rstart-w)).nonzero())
        nen     = np.max((rho<=(rfin+w)).nonzero())
        del n_used, rho, rstart, rfin, w_vals, w

        nreq   = 1 + (nen - nst)
        crange = np.arange(nreq)+nst
        self.rho_km_vals         = self.rho_km_vals[crange]
        self.start               = start-nst
        self.p_norm_vals         = self.p_norm_vals[crange]
        self.phase_rad_vals      = self.phase_rad_vals[crange]
        self.B_rad_vals          = self.B_rad_vals[crange]
        self.D_km_vals           = self.D_km_vals[crange]
        self.f_sky_hz_vals       = self.f_sky_hz_vals[crange]
        self.phi_rad_vals        = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals    = self.rho_dot_kms_vals[crange]
        self.T_hat_vals          = self.T_hat_vals[crange]
        self.F_km_vals           = self.F_km_vals[crange]
        self.w_km_vals           = self.w_km_vals[crange]
        self.mu_vals             = self.mu_vals[crange]
        self.lambda_sky_km_vals  = self.lambda_sky_km_vals[crange]
        del nreq, crange

    def __rect(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        w_func = np.zeros(nw_pts) + 1.0
        return w_func

    def __coss(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        w_func = np.cos(np.pi * x / w_in)**2
        return w_func

    def __kb20(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.0*np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb25(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5*np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb35(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 3.5 * np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kbmd20(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.0*np.pi
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    def __kbmd25(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.arange(nw_pts) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5 * np.pi
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func

    def __fresinv(self,T_hat,ker,dx,f_scale):
        T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
        return T

    def __fresinvfft(self,T_hat,ker,dx,f_scale):
        nw = np.size(T_hat)
        fft_t_hat       = np.fft.fft(T_hat)
        fft_conv        = np.fft.fft(ker)
        inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
        inv_t_hat      *= dx*(np.complex(1.0,1.0))/(2.0*f_scale)
        T               = inv_t_hat[int((nw-1)/2)]
        return T

    def __psifacfast(self,r,r0,cb,cp0,sp0):
        factor  = ((cb*cb) * cp0 * sp0 / (1.0 - (cb*cb) * (sp0*sp0))) * (r - r0) / r0
        return factor

    def __psifac(self,r,r0,b,phi0):
        cb      = np.cos(b)
        sp0     = np.sin(phi0)
        cp0     = np.cos(phi0)
        factor  = ((cb*cb) * cp0 * sp0 / (1.0 - (cb*cb) * (sp0*sp0))) * (r - r0) / r0
        return factor

    def __psifast(self,r,r0,d,cb,cp,sp,cp0,sp0):
        xi   = (cb / d) * (r0*cp0 - r*cp)
        eta  = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d*d)
        psi_vals   = np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi)
        return psi_vals

    def __normalize(self,r,w_func,f_scale):
        x         = r-np.mean(r)
        drho      = r[1]-r[0]
        f_scale   = f_scale
        psi       = (np.pi / 2.0) * ((x / f_scale)*(x / f_scale))
        ker       = np.exp(-1j * psi)
        T1        = np.abs(np.sum(w_func * ker) * drho)
        norm_fact = np.sqrt(2.0) * f_scale / T1
        return norm_fact

    def __fstnewrtraph(self,r,r0,d,b,phi,phi0):
        p1 = self.__psid1
        p2 = self.__psid2
        f1 = self.__psid1(r,r0,d,b,phi,phi0)
        f2 = self.__psid2(r,r0,d,b,phi,phi0)
        phi1 = phi - f1/f2
        f3 = self.__psid1(r,r0,d,b,phi1/f2,phi0)
        f4 = self.__psid2(r,r0,d,b,phi1/f2,phi0)
        y2 = phi-f1/f2-f3/f4
        return y2

    __func_dict = {
        "__rect" : {"func" : __rect, "normeq" : 1.00000000},
        "__coss" : {"func" : __coss, "normeq" : 1.50000000},
        "__kb20" : {"func" : __kb20, "normeq" : 1.49634231},
        "__kb25" : {"func" : __kb25, "normeq" : 1.65191895},
        "__kb35" : {"func" : __kb35, "normeq" : 1.92844639},
        "__kbmd20" : {"func" : __kbmd20, "normeq" : 1.52048174},
        "__kbmd25" : {"func" : __kbmd25, "normeq" : 1.65994218}
        }

    def __finv_n(self):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2. * np.pi * d_vals / lambda_vals
        cosb      = np.cos(B_rad_vals)
        cosphi0   = np.cos(phi_rad_vals)
        sinphi0   = np.sin(phi_rad_vals)
        dsq       = d_vals*d_vals
        rsq       = rho_vals*rho_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        psifac    = self.__psifacfast
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        psif      = self.__psifast
        nrm       = self.__normalize
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        w_init    = w_vals[start]
        w_func    = fw(w_init,dx)
        nw        = np.size(w_func)
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d      = d_vals[center]
            d2     = dsq[center]
            cb     = cosb[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            phi_s_rad1 = phi_s_rad
            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            
            # Compute psi and then compute the forward model.
            psi_vals = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
            F        = F_vals[center]
            # psi_vals = (np.pi/2.0)*(((r-r0)/F)*((r-r0)/F))
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center] = finv(T_hat,ker,dx,F)
            T_vals[center] *= nrm(r,w_func,F)
        return T_vals

    def __finv(self):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2. * np.pi * d_vals / lambda_vals
        cosb      = np.cos(B_rad_vals)
        cosphi0   = np.cos(phi_rad_vals)
        sinphi0   = np.sin(phi_rad_vals)
        dsq       = d_vals*d_vals
        rsq       = rho_vals*rho_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        psifac    = self.__psifacfast
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        psif      = self.__psifast
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        w_init    = w_vals[start]
        w_func    = fw(w_init,dx)
        nw        = np.size(w_func)
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d      = d_vals[center]
            d2     = dsq[center]
            cb     = cosb[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            phi_s_rad1 = phi_s_rad
            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            
            # Compute psi and then compute the forward model.
            psi_vals = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            F        = F_vals[center]
            T_vals[center] = finv(T_hat,ker,dx,F)
        return T_vals

    def __finv_n_v(self):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals     = 2. * np.pi * d_vals / lambda_vals
        cosb        = np.cos(B_rad_vals)
        cosphi0     = np.cos(phi_rad_vals)
        sinphi0     = np.sin(phi_rad_vals)
        dsq         = d_vals*d_vals
        rsq         = rho_vals*rho_vals
        # Define functions
        fw          = self.__func_dict[wtype]["func"]
        psifac      = self.__psifacfast
        if fft:
           finv     = self.__fresinvfft
        else: 
            finv    = self.__fresinv
        psif        = self.__psifast
        nrm         = self.__normalize
        # Calculate the corrected complex amplitude, point by point
        T_vals      = T_hat_vals * 0.0
        w_init      = w_vals[start]
        w_func      = fw(w_init,dx)
        nw          = np.size(w_func)
        pdb.set_trace()
        phi_s_rad1  = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d      = d_vals[center]
            d2     = dsq[center]
            cb     = cosb[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            phi_s_rad1 = phi_s_rad
            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            
            # Compute psi and then compute the forward model.
            psi_vals = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
            F        = F_vals[center]
            #psi_vals = (np.pi/2.0)*(((r-r0)/F)*((r-r0)/F))
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            
            T_vals[center] = finv(T_hat,ker,dx,F)
            T_vals[center] *= nrm(r,w_func,F)
            print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
            % (i,n_used,nw,loop),end="\r")
        return T_vals

    def __finv_v(self):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2.0 * np.pi * d_vals / lambda_vals
        cosb      = np.cos(B_rad_vals)
        cosphi0   = np.cos(phi_rad_vals)
        sinphi0   = np.sin(phi_rad_vals)
        dsq       = d_vals*d_vals
        rsq       = rho_vals*rho_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        psifac    = self.__psifacfast
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        psif      = self.__psifast
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        w_init    = w_vals[start]
        w_func    = fw(w_init,dx)
        nw        = np.size(w_func)
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d      = d_vals[center]
            d2     = dsq[center]
            cb     = cosb[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (np.abs(w_init - w)>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            phi_s_rad1 = phi_s_rad
            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            
            # Compute psi and then compute the forward model.
            psi_vals       = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
            ker            = w_func*np.exp(-1j*psi_vals)
            T_hat          = T_hat_vals[crange]
            F              = F_vals[center]
            T_vals[center] = finv(T_hat,ker,dx,F)
            print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
            % (i,n_used,nw,loop),end="\r")
        return T_vals

    def __finv_p(self,psitype):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2. * np.pi * d_vals / lambda_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            w_init = w
            w_func = fw(w,dx)
            nw     = np.size(w_func)
            crange = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r      = rho_vals[crange]
            F      = F_vals[center]
            x      = (r-r0)/F

            # Compute psi and then compute the forward model.
            psi_vals = (np.pi/2.0)*(x*x)
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            T_vals[center] = finv(T_hat,ker,dx,F)
        return T_vals

    def __finv_p_v(self,psitype):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2. * np.pi * d_vals / lambda_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            w_init = w
            w_func = fw(w,dx)
            nw     = np.size(w_func)
            crange = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r      = rho_vals[crange]
            F      = F_vals[center]
            x      = (r-r0)/F

            # Compute psi and then compute the forward model.
            psi_vals = (np.pi/2.0)*(x*x)
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            T_vals[center] = finv(T_hat,ker,dx,F)
            print("Pt: %d  Tot: %d  Width: %d Fast Inversion - Quadratic" \
            % (i,n_used,nw),end="\r")
        return T_vals

    def __finv_n_p(self,psitype):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2.0 * np.pi * d_vals / lambda_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        nrm       = self.__normalize
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            w_init = w
            w_func = fw(w,dx)
            nw     = np.size(w_func)
            crange = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r      = rho_vals[crange]
            F      = F_vals[center]
            B      = B_rad_vals[center]
            phi    = phi_rad_vals[center]
            D      = d_vals[center]
            kD     = kD_vals[center]
            afac   = np.cos(B)*np.cos(phi)
            if (psitype == 'taylor'):
                x      = (r-r0)/D
            else:
                x      = (r-r0)/F

            # Compute psi and then compute the forward model.
            if (psitype == 'taylor'):
                psi_vals    = (kD/2.0)*(1.0-afac*afac)*x*x
            else:
                psi_vals    = (np.pi/2.0)*x*x
            ker             = w_func*np.exp(-1j*psi_vals)
            T_hat           = T_hat_vals[crange]
            T_vals[center]  = finv(T_hat,ker,dx,F)
            #T_vals[center] *= nrm(r,w_func,F)
        return T_vals

    def __finv_n_p_v(self,psitype):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.D_km_vals
        B_rad_vals   = self.B_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        fft          = self.fft
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2. * np.pi * d_vals / lambda_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        nrm       = self.__normalize
        if fft:
           finv   = self.__fresinvfft
        else: 
            finv  = self.__fresinv
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            w_init = w
            w_func = fw(w,dx)
            nw     = np.size(w_func)
            crange = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r      = rho_vals[crange]
            F      = F_vals[center]
            x      = (r-r0)/F

            # Compute psi and then compute the forward model.
            psi_vals        = (np.pi/2.0)*(x*x)
            ker             = w_func*np.exp(-1j*psi_vals)
            T_hat           = T_hat_vals[crange]
            T_vals[center]  = finv(T_hat,ker,dx,F)
            T_vals[center] *= nrm(r,w_func,F)
            print("Pt: %d  Tot: %d  Width: %d Fast Inversion - Quadratic" \
            % (i,n_used,nw),end="\r")
        return T_vals

    def __write_tau_hist(self):
        """
        This creates a history dictionary.
        Returns:
            geo_hist (dict): Dictionary with "user name", "host name",
                    "run date", "python version", "operating system",
                    "source file", "input variables", and "input keywords".
        """
        user_name = os.getlogin()
        host_name = os.uname()[1]
        run_date  = time.ctime() + ' ' + time.tzname[0]
        python_v  = platform.python_version()
        opsys     = os.uname()[0]
        src_file  = __file__.split('/')[-1]
        src_dir   = __file__.rsplit('/',1)[0] +'/'
        rngreq    = (np.min(self.rng),np.max(self.rng))
        tau_hist = {
            "User Name"         : user_name,
            "Host Name"         : host_name,
            "Run Date"          : run_date,
            "Python Version"    : python_v,
            "Operating System"  : opsys,
            "Source Directory"  : src_dir,
            "Source File"       : src_file,
            "Input Variables"   : {
                "norm_inst": self.dathist,
                "res"      : self.res},
            "Input Keywords"    : {
                "b Factor"                      : self.bfac,
                "Resolution (km)"               : self.res,
                "Tapering Function"             : self.wtype,
                "Requested Radius (km)"         : rngreq,
                "Normalization"                 : self.norm,
                "Forward Model Data"            : self.fwd,
                "FFT's Used"                    : self.fft,
                "Psi Approximation"             : self.psitype},
        }
        return tau_hist

class compare_tau(object):
    def __init__(self,geodata,caldata,dlpdata,taudata,res,occ=False,
        rng='all',wtype="kb25",bfac=True,fft=False,verbose=True,
        norm=True):
        data    = extract_csv_data(geodata,caldata,dlpdata,
            occ=occ,taudata=taudata,verbose=verbose)
        
        tau_power   = data.power_vals
        tau_phase   = data.phase_vals
        tau_tau     = data.tau_vals
        rec         = diffraction_correction(data,res,rng=rng,bfac=bfac,
            wtype=wtype,fft=fft,verbose=verbose,norm=norm)
        tr          = data.tau_rho
        rho_km_vals = rec.rho_km_vals
        rmin        = np.min(tr)
        rmax        = np.max(tr)
        rstart      = int(np.min((rho_km_vals-rmin>=0).nonzero()))
        rfin        = int(np.max((rmax-rho_km_vals>=0).nonzero()))

        rho_km_vals = rec.rho_km_vals[rstart:rfin+1]
        power_vals  = rec.power_vals[rstart:rfin+1]
        phase_vals  = rec.phase_vals[rstart:rfin+1]
        tau_vals    = rec.tau_vals[rstart:rfin+1]

        rmin        = np.min(rho_km_vals)
        rmax        = np.max(rho_km_vals)
        rstart      = int(np.min((tr-rmin>=0).nonzero()))
        rfin        = int(np.max((rmax-tr>=0).nonzero()))

        tau_power   = tau_power[rstart:rfin+1]
        tau_phase   = tau_phase[rstart:rfin+1]
        tau_tau     = tau_tau[rstart:rfin+1]

        wtype       = rec.wtype
        res         = rec.res

        self.rho_km_vals = rho_km_vals
        self.power_vals  = power_vals
        self.tau_vals    = tau_vals
        self.phase_vals  = phase_vals
        self.tau_power   = tau_power
        self.tau_tau     = tau_tau
        self.tau_phase   = tau_phase
        self.wtype       = wtype
        self.res         = res

class find_optimal_resolution(object):
    def __init__(self,geo,cal,dlp,tau,occ,sres,dres,nres,
        rng='all',wlst=['kb25'],verbose=True):
        self.linfint = None
        self.l2int   = None
        self.linffft = None
        self.l2fft   = None

        nwins   = np.size(wlst)
        linfint = np.zeros((nwins,nres))
        l2int   = np.zeros((nwins,nres))
        linffft = np.zeros((nwins,nres))
        l2fft   = np.zeros((nwins,nres))
        resint  = np.zeros((nwins))
        resfft  = np.zeros((nwins))
        eres    = sres + (nres-1)*dres

        res         = sres
        data        = extract_csv_data(geo,cal,dlp,occ,taudata=tau,verbose=verbose)
        tr          = data.tau_rho
        tau_power   = data.power_vals
        rec         = diffraction_correction(data,res,rng=rng,wtype="kb35",verbose=False)
        rho_km_vals = rec.rho_km_vals
        power_vals  = rec.power_vals
        rmin        = np.min(rho_km_vals)
        rmax        = np.max(rho_km_vals)
        rstart      = int(np.min((tr-rmin>=0).nonzero()))
        rfin        = int(np.max((rmax-tr>=0).nonzero()))
        tau_power   = tau_power[rstart:rfin+1]
        sys.stdout.write("\033[K")
        for i in np.arange(nres):
            for j in range(nwins):
                wtype        = wlst[j]       
                recint       = diffraction_correction(data,res,rng=rng,wtype=wtype,fft=False)
                p_int        = recint.power_vals
                linf         = np.max(np.abs(p_int - tau_power))
                l2           = np.sqrt(np.sum(np.abs(p_int-tau_power)**2)*recint.dx_km)
                linfint[j,i] = linf
                l2int[j,i]   = l2

                recfft       = diffraction_correction(data,res,rng=rng,wtype=wtype,fft=True)
                p_fft        = recfft.power_vals
                linf         = np.max(np.abs(p_fft - tau_power))
                l2           = np.sqrt(np.sum(np.abs(p_fft-tau_power)**2)*recint.dx_km)
                linffft[j,i] = linf
                l2fft[j,i]   = l2
                sys.stdout.write("\033[K")
                if verbose:
                    printmes = ('Res:',res,'Max:',eres,"WTYPE:",wtype)
                    print("%s %f %s %f %s %s" % printmes)
            res += dres
        for j in range(nwins):
            resint[j] = sres+dres*np.min((linfint[j,...] == np.min(linfint[j,...])).nonzero())
            resfft[j] = sres+dres*np.min((linffft[j,...] == np.min(linffft[j,...])).nonzero())
        self.linfint = linfint
        self.l2int   = l2int
        self.linffft = linffft
        self.l2fft   = l2fft
        self.resint  = resint
        self.resfft  = resfft

class delta_impulse_diffraction(object):
    def __init__(self,geo,lambda_km,res,rho,dx_km_desired=0.25,
        occ=False,wtype='kb25',fwd=False,norm=True,bfac=True,
        verbose=True,fft=False,psitype='full',usefres=False):
    
        data = get_geo(geo,verbose=verbose)
        self.__retrieve_variables(data,verbose)
        self.__compute_variables(dx_km_desired,occ,verbose)
        self.__interpolate_variables(verbose)

        r0       = self.rho_km_vals
        b        = self.B_rad_vals
        d        = self.D_km_vals
        phi      = self.phi_rad_vals
        rng      = [rho-10.0*res,rho+10.0*res]
        nstar    = np.min((r0-rho>=0).nonzero())
        r        = r0[nstar]
        phistar  = phi[nstar]
        F        = fresnel_scale(lambda_km,d,phi,b)

        if not check_boole(fwd):
            raise TypeError("fwd must be Boolean: True/False")
        if not check_boole(norm):
            raise TypeError("norm must be Boolean: True/False")
        if not check_boole(bfac):
            raise TypeError("bfac must be Boolean: True/False")
        if not check_boole(fft):
            raise TypeError("fft must be Boolean: True/False")
        if not check_boole(verbose):
            raise TypeError("verbose must be Boolean: True/False")
        if not check_boole(usefres):
            raise TypeError("usefres must be Boolean: True/False")
        if (type(psitype) != type("Hi!")):
            raise TypeError("psitype must be a string: 'taylor', 'full', etc.")

        if (usefres == True):
            psi_vals = (np.pi/2.0)*((r0-r)/F)*((r0-r)/F)
        else:
            kD          = 2.0 * np.pi * d[nstar] / lambda_km
            dphi_s_rad  = psi_factor(r0,r,b,phi)
            phi_s_rad   = phi - dphi_s_rad
            cb          = np.cos(b)
            cp0         = np.cos(phi)
            sp0         = np.sin(phi)
            cp          = np.cos(phi_s_rad)
            sp          = np.sin(phi_s_rad)
            d2          = d*d
            r2          = r*r
            r02         = r0*r0

            loop        = 0
            dphi_s_rad1 = dphi_s_rad  
            dphi_s_rad  = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
                if verbose: print("Psi Iter: %d" % loop)
            psi_vals = kD * psi_fast(r,r0,d,cb,cp,sp,cp0,sp0)

        T_hat_vals     = (1.0-1.0j)*np.exp(1j*psi_vals)/(2.0*F)
        p_norm_vals    = np.abs(T_hat_vals)*np.abs(T_hat_vals)
        phase_rad_vals = -np.arctan2(np.imag(T_hat_vals),np.real(T_hat_vals))

        lambda_vals         = np.zeros(np.size(self.rho_km_vals))+lambda_km
        self.p_norm_vals    = p_norm_vals
        self.phase_rad_vals = phase_rad_vals
        self.f_sky_hz_vals  = freq_wav(lambda_vals)
        self.nstar          = nstar
        self.history        = "Bob"

        data = self

        recdata = diffraction_correction(data,res,rng=rng,
            wtype=wtype,fwd=fwd,norm=norm,verbose=verbose,
            bfac=bfac,fft=fft,psitype=psitype)

        self.rho_star           = rho
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.B_rad_vals         = recdata.B_rad_vals
        self.D_km_vals          = recdata.D_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rng                = recdata.rng
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq
        self.n_used             = recdata.n_used
        self.start              = recdata.start
        self.T_vals             = recdata.T_vals
        self.power_vals         = recdata.power_vals
        self.tau_vals           = recdata.tau_vals
        self.phase_vals         = recdata.phase_vals
        self.p_norm_fwd_vals    = None
        self.T_hat_fwd_vals     = None
        self.phase_fwd_vals     = None
        self.norm               = recdata.norm
        self.fwd                = recdata.fwd
        self.fft                = recdata.fft

    def __retrieve_variables(self,geo_dat,verbose):
        if verbose: print("Retrieving Variables...")
        geo_rho = np.array(geo_dat.rho_km_vals)
        rhotype = check_real(geo_rho)
        if not rhotype:
            raise TypeError("Bad GEO: rho_km_vals not real valued.")
        elif (np.min(geo_rho) < 0.0):
            raise ValueError("Bad GEO: rho_km_vals has negative values.")
        else: del rhotype

        geo_D   = np.array(geo_dat.D_km_vals)
        Dtype   = check_real(geo_D)
        if not Dtype:
            raise TypeError("Bad GEO: D_km_vals not real valued.")
        elif (np.min(geo_D) < 0.0):
            raise ValueError("Bad GEO: D_km_vals has negative values.")
        else: del Dtype

        geo_drho    = np.array(geo_dat.rho_dot_kms_vals)
        drhotype    = check_real(geo_drho)
        if not drhotype:
            raise TypeError("Bad GEO: rho_dot_kms_vals not real valued.")
        else: del drhotype

        geo_phi = np.array(geo_dat.phi_ora_deg_vals)
        phitype = check_real(geo_phi)
        if not phitype:
            raise TypeError("Bad GEO: phi_deg_ora_vals not real valued.")
        elif (np.max(np.abs(geo_phi)) > 360.0):
            raise ValueError("Bad GEO: phi_deg_ora_vals > 360")
        else: del phitype

        geo_B = np.array(geo_dat.B_deg_vals)
        Btype = check_real(geo_B)
        if not Btype:
            raise TypeError("Bad GEO: B_deg_vals not real valued.")
        elif (np.max(np.abs(geo_B)) > 360.0):
            raise ValueError("Bad GEO: B_de2g_vals > 360")
        else: del Btype

        self.geo_rho    = geo_rho
        self.geo_D      = geo_D
        self.geo_drho   = geo_drho
        self.geo_phi    = geo_phi
        self.geo_B      = geo_B

    def __compute_variables(self,dx_km_desired,occ,verbose):
        if verbose: print("Computing Variables...")
        geo_drho        = self.geo_drho

        if (occ == 'ingress'):  crange = (geo_drho < 0.0).nonzero()
        elif (occ == 'egress'): crange = (geo_drho > 0.0).nonzero()
        else:
            crange_e = (geo_drho > 0.0).nonzero()
            crange_i = (geo_drho < 0.0).nonzero()
            n_e      = np.size(crange_e)
            n_i      = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    "rho_dot_kms_vals has positive and negative values.\
                    This is likely a chord occultation. Set occ='ingress'\
                    to examine the ingress portion, and occ='egress'\
                    for the egress porition.")
            elif (n_e == 0) and (n_i == 0):
                raise ValueError("rho_dot_kms_vals is either empty or zero.")
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else: raise TypeError("Bad Input: GEO DATA")
            del n_e, n_i, crange_e, crange_i

        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else: raise ValueError("Bad occ input: Set 'egress' or 'ingress'")
            raise ValueError("Bad occ Input: '%s': %s" % (occ,mes))

        geo_rho             = self.geo_rho[crange]
        geo_D               = self.geo_D[crange]
        geo_drho            = self.geo_drho[crange]
        geo_phi             = self.geo_phi[crange]
        geo_B               = self.geo_B[crange]
        rmin                = np.min(geo_rho)
        rmax                = np.max(geo_rho)
        rho_km_vals         = np.arange(rmin,rmax,dx_km_desired)
        self.rho_km_vals    = rho_km_vals
        self.geo_rho        = geo_rho
        self.geo_D          = geo_D
        self.geo_drho       = geo_drho
        self.geo_phi        = geo_phi
        self.geo_B          = geo_B
        del rmin,rmax,geo_rho,geo_D,geo_drho,geo_phi,geo_B,rho_km_vals

    def __interpolate_variables(self,verbose):
        if verbose: print("Interpolating Data...")
        rho_km_vals      = self.rho_km_vals
        geo_rho          = self.geo_rho
        geo_drho         = self.geo_drho
        geo_D            = self.geo_D
        geo_phi          = self.geo_phi
        geo_B            = self.geo_B
        D_km_interp      = interpolate.interp1d(geo_rho,geo_D,kind='linear')
        D_km_vals        = D_km_interp(rho_km_vals)
        rho_dot_interp   = interpolate.interp1d(geo_rho,geo_drho,kind='linear')
        rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        phi_deg_interp   = interpolate.interp1d(geo_rho,geo_phi,kind='linear')
        phi_deg_vals     = phi_deg_interp(rho_km_vals)
        B_deg_interp     = interpolate.interp1d(geo_rho,geo_B,kind='linear')
        B_deg_vals       = B_deg_interp(rho_km_vals)
        phi_rad_vals     = np.deg2rad(phi_deg_vals)
        B_rad_vals       = np.deg2rad(B_deg_vals)

        del geo_rho,geo_drho,geo_D,geo_phi,geo_B,D_km_interp,rho_dot_interp
        del phi_deg_vals,phi_deg_interp,B_deg_interp,B_deg_vals
        
        self.rho_km_vals      = rho_km_vals
        self.phi_rad_vals     = phi_rad_vals
        self.B_rad_vals       = B_rad_vals
        self.D_km_vals        = D_km_vals
        self.rho_dot_kms_vals = rho_dot_kms_vals