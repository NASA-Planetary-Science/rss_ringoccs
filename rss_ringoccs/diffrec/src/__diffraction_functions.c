/******************************************************************************
 *                          Diffraction Functions                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains routines to compute the Fresnel transform for      *
 *      diffraction reconstruction or forward modeling.                       *
 *          Fresnel Quadratic Approximations:                                 *
 *              Classic quadratic approximation used in Fourier optics.       *
 *          Legendre Expansion:                                               *
 *             Approximate the Fresnel kernel using Legendre polynomials.     *
 *          Newton-Raphson Method:                                            *
 *              Performs the Fresnel Transform by computing the stationary    *
 *              value of the Fresnel Kernel using the Newton-Raphson method.  *
 ******************************************************************************
 *  The Inverse Fresnel Transform:                                            *
 *                                                                            *
 *                infinity                                                    *
 *                     -                                                      *
 *                    | |                                                     *
 *         T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0           *
 *                  | |                                                       *
 *                   -                                                        *
 *               -infinity                                                    *
 *                                                                            *
 *      Where T_hat is the diffracted data, w is the window function, r is    *
 *      the ring intercept point, and r_0 is a dummy variable of integration. *
 *      psi is the Fresnel Kernel, and exp is the exponential function.       *
 ******************************************************************************
 *  The Normalization Scheme:                                                 *
 *      As the resolution get's too high, say 10 km or larger, the window     *
 *      width quickly shrinks to zero and the integral will be approximately  *
 *      zero. To account for this, the option to normalize the integral by    *
 *      the window width is offered. The normalization is defined as follows: *
 *                                                                            *
 *                    |     _ +infinity           |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -  -infinity            |                           *
 *          Norm =  __________________________________                        *
 *                  |    -  +W/2                    |                         *
 *                  |   | |                         |                         *
 *                  |   |    w(x) exp(-i psi(x)) dx |                         *
 *                  | | |                           |                         *
 *                  |  -   -W/2                     |                         *
 *                                                                            *
 *      This has the effect of making the free-space regions, which are       *
 *      regions that are not affected by diffraction, evaluate to one,        *
 *      regardless of what (positive) resolution was chosen.                  *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  reset_window:                                                             *
 *      Void function that takes in a pointer to a double array and creates   *
 *      an array of values for the ring radius within half a window width of  *
 *      the ring intercept point. Do to symmetry, only the values to the left *
 *      of the ring intercept point are computed. This also computes the      *
 *      window function as a function of this array.                          *
 ******************************************************************************
 *  Fresnel_Transform_Double / Fresnel_Transform_Norm_Double                  *
 *      Computes the Fresnel transform using a quadratic approximation to the *
 *      Fresnel kernel, with or without normalization.                        *
 ******************************************************************************
 *  Fresnel_Legendre_Double / Fresnel_Legendre_Norm_Double                    *
 *      Computes the Fresnel Transform using Legendre Polynomials to          *
 *      approximate the Fresnel kernel to various powers (4, 6, or 8).        *
 ******************************************************************************
 *  Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double    *
 *      Computes the Fresnel inverse transform using Newton-Raphson to        *
 *      compute the stationary value of the Fresnel kernel.                   *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code uses complex numbers throughout, and is compatible with the C99 *
 *  standard. To use this code, make sure your compiler supports C99 or more  *
 *  recent standards of the C Programming Language.                           *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************/
#include "__diffraction_functions.h"

/******************************************************************************
 *------------------------------DEFINE C FUNCTIONS----------------------------*
 *  These are functions written in pure C without the use of the C-Python API.*
 ******************************************************************************/

/******************************************************************************
 *  Function:                                                                 *
 *      reset_window                                                          *
 *  Purpose:                                                                  *
 *      Compute an array of points from -width to zero, equally spaced by dx. *
 *      This acts as the independent variable for later use.                  *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed, and rho is the dummy variable of integration *
 *          which varies from rho0-W/2 to rho0, W being the window width.     *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      width (double):                                                       *
 *          The width of the window function.                                 *
 *      nw_pts (long):                                                        *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      fw  double (*)(double, double):                                       *
 *          Function pointer to the window function.                          *
 *  Notes:                                                                    *
 *      1.) This is a void function that takes in pointers as arguments. The  *
 *          values of the pointers are changed within this function and there *
 *          is no need to return anything. Hence, no return statement.        *
 ******************************************************************************/
static void reset_window(double *x_arr, double *w_func, double dx, double width,
                         long nw_pts, double (*fw)(double, double))
{
    /*  Create a variable for indexing.                                       */
    long i;

    /* Loop over i, computing the window function and the x_arr variable.     */
    for(i=0; i<nw_pts; ++i){
        x_arr[i] = (i-nw_pts)*dx;
        w_func[i] = fw(x_arr[i], width);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Double / Fresnel_Transform_Norm_Double              *
 *  Purpose:                                                                  *
 *      Perform diffraction correction on diffraction limited data using the  *
 *      classic Fresnel quadratic approximation to the Fresnel kernel.        *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as pi/2 (rho-rho0)^2, where rho is the ring radius of     *
 *          the point being reconstructed, and rho0 is the dummy variable of  *
 *          integration which varies from rho-W/2 to zero, W being the        *
 *          window width. This should point to n_pts number of elements.      *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr. This should *
 *          point to n_pts number of elements.                                *
 *      F (double):                                                           *
 *          The Fresnel scale at rho.                                         *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
complex double Fresnel_Transform_Double(double *x_arr, complex double *T_in,
                                        double *w_func, double F, double dx,
                                        long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2;

    /*  exp_negative_ix is used for the Fresnel kernel.                       */
    complex double T_out, exp_negative_ix;

    /*  Initialize the T_out variable to zero, so we can loop over later.     */
    T_out = 0.0;

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    j = n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprocal of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/F;
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        T_out += exp_negative_ix * (T_in[center - j] + T_in[center + j]);
        j -= 1;
    }

    /*  Add the central point in the Riemann sum. This is the center of the   *
     *  window function. That is, where w_func = 1.                           */
    T_out += T_in[center];

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5+0.5*_Complex_I)*dx*rcpr_F;
    return T_out;
}

complex double Fresnel_Transform_Norm_Double(double *x_arr,
                                             complex double *T_in,
                                             double *w_func, double F,
                                             double dx, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2;

    /*  exp_negative_ix is the Fresnel kernel, norm is the normalization.     */
    complex double T_out, exp_negative_ix, norm;

    /*  Initialize T_out and norm to zero, so we can loop over later.         */
    T_out  = 0.0;
    norm   = 0.0;

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    j = n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprical of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/F;
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += 2.0*exp_negative_ix;

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        T_out += exp_negative_ix * (T_in[center + j] + T_in[center - j]);
        j -= 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    T_out += T_in[center];
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5+0.5*_Complex_I)*norm;
    return T_out;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Legendre_Double / Fresnel_Legendre_Double_Norm                *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using Legendre polynomials to         *
 *      approximate the Fresnel kernel. Do to the nature of Legendre          *
 *      polynomials and the first iteration of the Newton-Raphson scheme      *
 *      applied to the Fresnel kernel, this is often extremely accurate and   *
 *      extremely fast.                                                       *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed and rho is the dummy variable of integration  *
 *          which varies from rho0-W/2 to rho0, W being the window width.     *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr. This should *
 *          contain n_pts (see below) number of elements.                     *
 *      D (double):                                                           *
 *          The distance from the spacecraft to the ring-intercept point.     *
 *      coeffs (double *):                                                    *
 *          The coefficients of the polynomial approximation of psi. This     *
 *          should conter order (see below) number of elements.               *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      F (double):                                                           *
 *          The Fresnel scale, in the same units as D and dx.                 *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      order (int):                                                          *
 *          Defined as degree-1, where degree is the degree of the polynomial *
 *          approximation for psi. Order is the highest Legendre polynomial   *
 *          to be used. This is also the size of coeffs (See above).          *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
complex double Fresnel_Legendre_Double(double *x_arr, complex double *T_in,
                                       double *w_func, double D, double *coeffs,
                                       double dx, double F, double kd,
                                       long n_pts, unsigned char order,
                                       long center)
{
    /*  Declare all necessary variables. i, j, and k are used for indexing.   */
    long i, j, k;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi;
    double psi_even, psi_odd;
    double rcpr_D;
    complex double exp_negative_psi, exp_positive_psi;
    complex double T_out;

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    rcpr_D = 1.0/D;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x  = x_arr[i]*rcpr_D;
        x2 = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = coeffs[order-1];
        psi_odd = coeffs[order-2];
        for (k=3; k<order-1; ++k){
            psi_even = psi_even*x2 + coeffs[order-k];
            psi_odd  = psi_odd*x2 + coeffs[order-k-1];
        }

        /*  The leading term is x^2, so multiply by this and kD.              */
        psi_even = psi_even*x2 + coeffs[0];
        psi_even *= kd*x2;
        psi_odd *= kd*x2*x;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        T_out += exp_negative_psi * T_in[center - j];
        T_out += exp_positive_psi * T_in[center + j];
        j -= 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    T_out += T_in[center];

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx / F;
    return T_out;
}

complex double Fresnel_Legendre_Norm_Double(double *x_arr, complex double *T_in,
                                            double *w_func, double D,
                                            double *coeffs, double dx, double F,
                                            double kd, long n_pts,
                                            unsigned char order, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j, k;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi;
    double psi_even, psi_odd;
    double rcpr_D;
    complex double exp_negative_psi, exp_positive_psi;
    complex double T_out;
    complex double norm;

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    rcpr_D = 1.0/D;


    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;
    norm  = 0.0;

    j = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x  = x_arr[i]*rcpr_D;
        x2 = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.  */
        psi_even = coeffs[order-1];
        psi_odd = coeffs[order-2];
        for (k=3; k<order-1; ++k){
            psi_even = psi_even*x2 + coeffs[order-k];
            psi_odd  = psi_odd*x2 + coeffs[order-k-1];
        }

        /*  The leading term is x^2, so multiply by this and kD.              */
        psi_even  = psi_even*x2 + coeffs[0];
        psi_even *= kd*x2;
        psi_odd  *= kd*x2*x;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += exp_negative_psi+exp_positive_psi;

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        T_out += exp_negative_psi * T_in[center - j];
        T_out += exp_positive_psi * T_in[center + j];
        j -= 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    T_out += T_in[center];
    norm += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double*
 *  Purpose:                                                                  *
 *      Performs diffraction correction by using the Newton-Raphson method of *
 *      root-finding to compute the stationary value of the Fresnel kernel.   *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed and rho is the dummy variable of integration  *
 *          which varies from rho0-W/2 to rho+W/2, W being the window width.  *
 *          This should contain n_pts (see below) number of elements.         *
 *      phi_arr (double *):                                                   *
 *          The ring azimuth angle corresponding to rho, in radians.          *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *          This should contain n_pts (see below) number of elements.         *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      r (double):                                                           *
 *          The radius of the current point be reconstructed.                 *
 *      B (double):                                                           *
 *          The ring opening angle, in radians.                               *
 *      D (double):                                                           *
 *          The distance from the spacecraft to the ring-intercept point.     *
 *      EPS (double):                                                         *
 *          The allowable error in the computation of the stationary value of *
 *          the Fresnel kernel.                                               *
 *      toler (long):                                                         *
 *          The maximum number of iterations the Newton-Raphson method is     *
 *          allowed to undergo.                                               *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      F (double):                                                           *
 *          The Fresnel scale, in the same units as D and dx.                 *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
complex double Fresnel_Transform_Newton_Double(
    double *x_arr, double *phi_arr, complex double *T_in, double *w_func,
    double kD, double r, double B, double D, double EPS, long toler, double dx,
    double F, long n_pts, long center
){

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi;
    complex double T_out, exp_psi;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi     = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx / F;
    return T_out;
}

complex double Fresnel_Transform_Newton_Norm_Double(
    double *x_arr, double *phi_arr, complex double *T_in, double *w_func,
    double kD, double r, double B, double D, double EPS, long toler, double dx,
    double F, long n_pts, long center
){

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi;
    complex double T_out, exp_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi     = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the norm using a Riemann sum as well.                     */
        norm   += exp_psi;

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

void DiffractionCorrectionFresnel(DLPObj dlp)
{
    /*  i and j used for indexing, nw_pts is number of points in window.     */
    long i, j, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).  */
    double w_init, dx, two_dx;

    /*  Declare function pointers. fw is window function, FresT is Fresnel   *
     *  transform. Options for FresT are with or without normalization.      */
    double (*fw)(double, double);
    complex double (*FresT)(double*, complex double*, double*, double,
                            double, long, long);

    /*  Cast the selected window type to the fw pointer.                     */
    if      (dlp.wtype == 0){fw = &Rect_Window_Double;}
    else if (dlp.wtype == 1){fw = &Coss_Window_Double;}
    else if (dlp.wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (dlp.wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (dlp.wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else                    {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    /*  Cast FresT to the appropriate function.                               */
    if (dlp.use_norm){FresT = &Fresnel_Transform_Norm_Double;}
    else {FresT = &Fresnel_Transform_Double;}

    /* Move the pointers to the correct starting point. Compute window width.*/
    center  = dlp.start;

    /*  Compute some extra necessary variables.                               */
    w_init  = dlp.w_km_vals[center];
    dx      = dlp.rho_km_vals[center+1] - dlp.rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx.                                          */
    double* x_arr  = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);

    /*  Pass the x_arr array (ring radius) to the void function get_arr.      *
     *  This alters the x_arr pointer so that it's values range from -W/2 to  *
     *  zero, where W is the window width.                                    */
    reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Compute Window Functions, and compute pi/2 * x^2                       */
    for(j=0; j<nw_pts; ++j){
        x_arr[j] *= PI_BY_TWO*x_arr[j];

        /*  If forward transform is selected, negate x_arr.                   */
        x_arr[j] *= (dlp.use_fwd == 0) - (dlp.use_fwd == 1);
    }

    /*  Compute the Fresnel transform across the input data.                  */
    for (i=0; i<=dlp.n_used; ++i){
        if (fabs(w_init - dlp.w_km_vals[center]) >= two_dx) {
            /* Reset w_init and recompute window function.                    */
            w_init = dlp.w_km_vals[center];
            nw_pts = (long)(w_init / two_dx);

            /*  Reallocate memory, since the sizes of the arrays changed.     */
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.         */
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

            /* Compute Window Functions, and compute pi/2 * x^2               */
            for(j=0; j<nw_pts; ++j){
                x_arr[j] *= PI_BY_TWO*x_arr[j];
                x_arr[j] *= (dlp.use_fwd == 0) - (dlp.use_fwd == 1);
            }
        }

        /*  Compute the Fresnel Transform about the current point.           */
        dlp.T_out[center] = FresT(x_arr, dlp.T_in, w_func,
                                  dlp.F_km_vals[center], dx, nw_pts, center);

        /*  Move the pointers to the next point.                             */
        center += 1;
    }
    free(x_arr);
    free(w_func);
}

void DiffractionCorrectionLegendre(DLPObj dlp)
{
    long i, nw_pts, center;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff;
    double (*fw)(double, double);
    complex double (*FresT)(double*, complex double*, double*, double, double *,
                            double, double, double, long, unsigned char, long);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp.wtype == 0){fw = &Rect_Window_Double;}
    else if (dlp.wtype == 1){fw = &Coss_Window_Double;}
    else if (dlp.wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (dlp.wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (dlp.wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else                    {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    if (dlp.use_norm){FresT = &Fresnel_Legendre_Norm_Double;}
    else             {FresT = &Fresnel_Legendre_Double;}

    /* Compute first window width and window function. */
    center = dlp.start;

    if (dlp.use_fwd){
        for (i=0; i <= dlp.n_used; ++i){
            dlp.kd_vals[i] *= -1.0;
        }
    }

    w_init  = dlp.w_km_vals[center];
    dx      = dlp.rho_km_vals[center+1] - dlp.rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    double* x_arr              = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func             = (double *)malloc(sizeof(double)*nw_pts);
    double* legendre_p         = (double *)malloc(sizeof(double)*(dlp.order+1));
    double* alt_legendre_p     = (double *)malloc(sizeof(double)*dlp.order);
    double* fresnel_ker_coeffs = (double *)malloc(sizeof(double)*dlp.order);

    reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    for (i = 0; i <= dlp.n_used; ++i){
        cosb            = cos(dlp.B_rad_vals[center]);
        cosp            = cos(dlp.phi_rad_vals[center]);
        sinp            = sin(dlp.phi_rad_vals[center]);
        Legendre_Coeff  = cosb*sinp;
        Legendre_Coeff *= Legendre_Coeff;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        /* Compute Legendre Polynomials,                                      */
        Legendre_Polynomials(legendre_p, cosb*cosp, dlp.order+1);
        Alt_Legendre_Polynomials(alt_legendre_p, legendre_p, dlp.order);

        /*  Compute the coefficients using Cauchy Products. First compute     *
         *  the bottom triangle of the square in the product.                 */
        Fresnel_Kernel_Coefficients(fresnel_ker_coeffs, legendre_p,
                                    alt_legendre_p, Legendre_Coeff, dlp.order);

        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - dlp.w_km_vals[center]) >= two_dx) {

            /* Reset w_init and recompute window function.                    */
            w_init = dlp.w_km_vals[center];
            nw_pts = (long)(w_init / two_dx);
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);
        }

        /*  Compute the fresnel tranform about the current point.             */
        dlp.T_out[center] = FresT(x_arr, dlp.T_in, w_func,
                                  dlp.D_km_vals[center], fresnel_ker_coeffs, dx,
                                  dlp.F_km_vals[center], dlp.kd_vals[center],
                                  nw_pts, dlp.order, center);

        /*  Increment T_in pointer using pointer arithmetic.                  */
        center += 1;
    }

    free(x_arr);
    free(w_func);
    free(legendre_p);
    free(alt_legendre_p);
    free(fresnel_ker_coeffs);
}

void DiffractionCorrectionNewton(DLPObj dlp)
{
    long i, j, nw_pts, toler, center;
    double w_init, dx, two_dx, EPS;

    toler = 5;
    EPS = 1.E-4;

    double (*fw)(double, double);
    complex double (*FresT)(double *, double *, complex double *, double *,
                            double, double, double, double, double, long,
                            double, double, long, long);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp.wtype == 0){fw = &Rect_Window_Double;}
    else if (dlp.wtype == 1){fw = &Coss_Window_Double;}
    else if (dlp.wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (dlp.wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (dlp.wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (dlp.wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else                     {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    if (dlp.use_norm){FresT = &Fresnel_Transform_Newton_Norm_Double;}
    else {FresT = &Fresnel_Transform_Newton_Double;}

    /* Compute first window width and window function. */
    center = dlp.start;

    if (dlp.use_fwd){
        for (i=0; i<=dlp.n_used; ++i){
            dlp.kd_vals[center+i] *= -1.0;
        }
    }

    w_init  = dlp.w_km_vals[center];
    dx      = dlp.rho_km_vals[center+1] - dlp.rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    double *x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double *phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double *w_func  = (double *)malloc(sizeof(double) * nw_pts);

    for (j=0; j<nw_pts; ++j){
        x_arr[j]   = dlp.rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = dlp.phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - dlp.rho_km_vals[center], w_init);
    }

    for (i=0; i<=dlp.n_used; ++i){

        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - dlp.w_km_vals[center]) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = dlp.w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = dlp.rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp.phi_rad_vals[center+j-(nw_pts-1)/2];
                w_func[j]  = fw(x_arr[j] - dlp.rho_km_vals[center], w_init);
            }
        }
        else {
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = dlp.rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp.phi_rad_vals[center+j-(nw_pts-1)/2];
            }
        }

        /*  Compute the fresnel tranform about the current point.   */
        dlp.T_out[center] = FresT(x_arr, phi_arr, dlp.T_in, w_func,
                                  dlp.kd_vals[center], dlp.rho_km_vals[center],
                                  dlp.B_rad_vals[center], dlp.D_km_vals[center],
                                  EPS, toler, dx, dlp.F_km_vals[center],
                                  nw_pts, center);

        /*  Increment pointers using pointer arithmetic.                      */
        center += 1;
    }
    free(x_arr);
    free(phi_arr);
    free(w_func);
}