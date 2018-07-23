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