/*******************************************************************************
 *                          Diffraction Functions                              *
 *******************************************************************************
 *  Purpose:                                                                   *
 *      This file contains functions used for computing the Fresnel Inverse    *
 *      Transform on a set of diffraction limited data. There are several      *
 *      Methods of performing this:                                            *
 *          Fresnel Quadratic Approximations:                                  *
 *              Classic quadratic approximation that is used in Fourier Optics.*
 *          Legendre Cubic Expansion:                                          *
 *              Cubic approimation of Fresnel Kernel by Legenedre Polynomials. *
 *          Legendre Quartic Expansion:                                        *
 *              Quartic approximation of Fresnel Kernel                        *
 *              by Legendre polynomials.                                       *
 *          Legendre Sextic Expansion:                                         *
 *              Sextic approximation of Fresnel Kernel by Legendre Polynomials.*
 *          Legendre Octic Expansion:                                          *
 *              Octic approximation of Fresnel Kernel by Legendre Polynomials. *
 *          Newton-Raphson Method:                                             *
 *              Computes the Fresnel Inverse Transform by computing the        *
 *              stationary value of the Fresnel Kernel using the               *
 *              Newton-Raphson method of root-finding.                         *
 *******************************************************************************
 *  Variables:                                                                 *
 *     A_0         (Double):                                                   *
 *         The coefficient of the x^2 term in the expansion for psi.           *
 *     A_1         (Double):                                                   *
 *         The coefficient of the x^3 term in the expansion for psi.           *
 *     A_2         (Double):                                                   *
 *         The coefficient of the x^4 term in the expansion for psi.           *
 *     A_3         (Double):                                                   *
 *         The coefficient of the x^5 term in the expansion for psi.           *
 *     A_4         (Double):                                                   *
 *         The coefficient of the x^6 term in the expansion for psi.           *
 *     A_5         (Double):                                                   *
 *         The coefficient of the x^7 term in the expansion for psi.           *
 *     A_6         (Double):                                                   *
 *         The coefficient of the x^8 term in the expansion for psi.           *
 *     dx          (Double):                                                   *
 *         Spacing between two consecutive ring intercept points, km.          *
 *     kd          (Double):                                                   *
 *         The wavenumber, k, weighted by the spacecraft-ring distance, D.     *
 *     n_pts       (Long):                                                     *
 *         Half the number of points in the window, rounded down.              *
 *     rcpr_D      (Double):                                                   *
 *         1/D, where D is the distant from the spacecraft to the ring         *
 *         intercept point, in kilometers.                                     *
 *     rcpr_F      (Double):                                                   *
 *         The reciprocal of the Fresnel scale in kilometers.                  *
 *     T_in        (Pointer to Char):                                          *
 *         The raw diffraction-limited data that is to be corrected.           *
 *     T_in_steps  (npy_intp (Equivalent to Long)):                            *
 *         The number of steps in memory to get from the nth data point        *
 *         to the (n+1)th data point in the T_in variable (See above).         *
 *     T_out       (Complex Double):                                           *
 *         The diffraction corrected profile.                                  *
 *     w_func      (Pointer to Double):                                        *
 *         Pre-computed window function. Should have the same number of        *
 *         points as the x_arr pointer.                                        *
 *     x_arr       (Pointer to Double):                                        *
 *         The ring radii within the given window.                             *
 *******************************************************************************
 *  The Inverse Fresnel Transform:                                             *
 *                                                                             *
 *                infinity                                                     *
 *                     -                                                       *
 *                    | |                                                      *
 *         T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0            *
 *                  | |                                                        *
 *                   -                                                         *
 *               -infinity                                                     *
 *                                                                             *
 *      Where T_hat is the diffracted data, w is the window function, r is     *
 *      the ring intercept point, and r_0 is a dummy variable of integration.  *
 *      psi is the Fresnel Kernel, and exp is simply the exponential function. *
 *******************************************************************************
 *  The Normalization Scheme:                                                  *
 *      As the resolution get's too high, say 10 km or greater, the window     *
 *      width quickly shrinks to zero. Thus the integral will be approximately *
 *      zero. To account for this, the option to normalize the integral by the *
 *      window width is offered. The normalization is defined as follows:      *
 *                                                                             *
 *                    |     _ +infinity           |                            *
 *                    |    | |                    |                            *
 *                    |    |    exp(-i psi(x)) dx |                            *
 *                    |  | |                      |                            *
 *                    |   -  -infinity            |                            *
 *          Norm =  __________________________________                         *
 *                  |    -  +W/2                    |                          *
 *                  |   | |                         |                          *
 *                  |   |    w(x) exp(-i psi(x)) dx |                          *
 *                  | | |                           |                          *
 *                  |  -   -W/2                     |                          *
 *                                                                             *
 *      This has the effect of making the free-space regions, or regions which *
 *      were not affected by diffraction, evaluate to approximately one,       *
 *      regardless of what resolution was chosen.                              *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  get_arr:                                                                   *
 *      Void function that takes in a pointer to a double array and creates    *
 *      an array of values for the ring radius within half a window width of   *
 *      the ring intercept point. Do to symmetry, only the values to the left  *
 *      of the ring intercept point are computed.                              *
 *******************************************************************************
 *  _fresnel_transform:                                                        *
 *      Computes the Fresnel Inverse Transform using the classic Fresnel       *
 *      quadratic approximation. No normalization is applied.                  *
 *******************************************************************************
 *  _fresnel_transform_norm:                                                   *
 *      Same as _fresnel_transform, but the normalization is applied.          *
 *******************************************************************************
 *  _fresnel_cubic,                                                            *
 *  _fresnel_quartic,                                                          *
 *  _fresnel_sextic,                                                           *
 *  _fresnel_octic:                                                            *
 *      Computes the Fresnel Inverse Transform using Legendre Polynomials to   *
 *      approximate the Fresnel kernel to various powers (3, 4, 6, or 8).      *
 *******************************************************************************
 *  _fresnel_cubic_norm,                                                       *
 *  _fresnel_quartic_norm,                                                     *
 *  _fresnel_sextic_norm,                                                      *
 *  _fresnel_octic_norm:                                                       *
 *      Same as previous functions, but with the normalization scheme.         *
 *******************************************************************************
 *  _fresnel_transform_newton:                                                 *
 *      Computes the Fresnel inverse transform using Newton-Raphson to compute *
 *      the stationary value of the Fresnel kernel.                            *
 *******************************************************************************
 *  fresnel_transform_newton_norm:                                             *
 *      Same as previous function, but with the normalization scheme.          *
 *******************************************************************************
 *                            A FRIENDLY WARNING                               *
 *******************************************************************************
 *  This code uses complex numbers throughout, and is compatible with the C99  *
 *  standard. To use this code, make sure your compiler supports C99 or more   *
 *  recent standards of the C Programming Language.                            *
 *******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                                *
 *  Date:       June 21, 2019                                                  *
 ******************************************************************************/

/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H
#define RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H

/*  Various trig functions, complex variables, and more found here.           */
#include <math.h>
#include <complex.h>

/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  Functions for computing the Fresnel Kernel and Newton's Method.           */
#include "__fresnel_kernel.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 *  These are functions written in pure C without the use of the Numpy-C API.  *
 *  They are used to define various special functions. They will be wrapped in *
 *  a form that is useable with the Python interpreter later on within         *
 *  the _diffraction_functions.c file.                                         *
 ******************************************************************************/

/*  This void function takes a pointer as an input. The pointers values       *
 *  are changed within the code, and there is no need to "return" anything.   *
 *  Hence, this function has no "return" statement.                           */
static void get_arr(double* x_arr, double dx, long nw_pts)
{
    /***************************************************************************
     *  Function:                                                              *
     *      get_arr                                                            *
     *  Purpose:                                                               *
     *      This computes an array of length nw_pts, with values ranging from  *
     *      -nw_pts*dx to zero. Due to symmetry in the reconstruction, only    *
     *      the left half of a given window is needed, and hence this returns  *
     *      only half of the array. The values zero to nw_pts*dx are omitted.  *
     **************************************************************************/

    /*  Declare variable for indexing over.                                   */
    long i;

    /*  Loop over the input array, and assign the correct ring radius (km)    */
    for (i=0; i<nw_pts; ++i){
        x_arr[i] = (i-nw_pts)*dx;
    }
}

complex double _fresnel_transform(double* x_arr, char* T_in, double* w_func,
                                  double F, double dx, long n_pts,
                                  npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform                                                 *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses the classic Fresnel quadratic aproximation to   *
     *      the Fresnel Kernel to perform diffraction correction for the given *
     *      input data. Mathematical definitions are given in the comment      *
     *      at the start of this file.                                         *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2;

    /*  exp_negative_ix is used for the Fresnel kernel.                       */
    complex double T_out, exp_negative_ix;

    /*  Initialize the T_out variable to zero, so we can loop over later.     */
    T_out = 0.0;

    j = -n_pts;

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
         *  This cuts the number of computations roughly in half.             */
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is the center of the   *
     *  window function. That is, where w_func = 1.                           */
    T_out += *(complex double *)T_in;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5+0.5*_Complex_I)*dx*rcpr_F;
    return T_out;
}

complex double _fresnel_transform_norm(double* x_arr, char* T_in,
                                       double* w_func, double F, double dx,
                                       long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_norm                                            *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses the classic Fresnel quadratic aproximation to   *
     *      the Fresnel Kernel to perform diffraction correction for the given *
     *      input data, and applies the normalization scheme to the output.    *
     *      Mathematical definitions are given in the comment at the start of  *
     *      this file.                                                         *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2;

    /*  exp_negative_ix is the Fresnel kernel, norm is the normalization.     */
    complex double T_out, exp_negative_ix, norm;

    /*  Initialize T_out and norm to zero, so we can loop over later. */
    T_out  = 0.0;
    norm   = 0.0;

    j = -n_pts;

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
         *  This cuts the number of computations roughly in half.             */
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5+0.5*_Complex_I)*norm;
    return T_out;
}

complex double _fresnel_cubic(double* x_arr, char* T_in, double* w_func,
                              double rcpr_D, double A_0, double A_1, double dx,
                              double rcpr_F, double kd, long n_pts,
                              npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_cubic                                           *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a cubic expansion, and then performs diffraction      *
     *      correction on diffracted data using this. Mathematical definitions *
     *      are given in the comment at the start of this file.                *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*A_0;
        psi_odd  = kd*x2*x*A_1;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.                      */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

complex double _fresnel_cubic_norm(double* x_arr, char* T_in, double* w_func,
                                   double rcpr_D, double A_0, double A_1,
                                   double dx, double rcpr_F, double kd,
                                   long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_cubic_norm                                      *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a cubic expansion, and then performs diffraction      *
     *      correction on diffracted data using this. The normalization scheme *
     *      is applied to the output. Mathematical definitions are given in    *
     *      the comment at the start of this file.                             *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*A_0;
        psi_odd  = kd*x2*x*A_1;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += exp_negative_psi+exp_positive_psi;

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double _fresnel_quartic(double* x_arr, char* T_in, double* w_func,
                                double rcpr_D, double A_0, double A_1,
                                double A_2, double dx, double rcpr_F, double kd,
                                long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_quartic                                         *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this. Mathematical definitions *
     *      are given in the comment at the start of this file.                *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*A_2);
        psi_odd  = kd*x2*x*A_1;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

complex double _fresnel_quartic_norm(double* x_arr, char* T_in, double* w_func,
                                     double rcpr_D, double A_0, double A_1,
                                     double A_2, double dx, double rcpr_F,
                                     double kd, long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_quartic                                         *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this. The normalization scheme *
     *      is then applied to the output. Mathematical definitions can be     *
     *      found in the comment at the start of this file.                    *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*A_2);
        psi_odd  = kd*x2*x*A_1;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += exp_negative_psi+exp_positive_psi;

        /*  Approximate the integral with a Riemann Sum.                      */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double _fresnel_sextic(double* x_arr, char* T_in, double* w_func,
                               double rcpr_D, double A_0, double A_1,
                               double A_2, double A_3, double A_4, double dx,
                               double rcpr_F, double kd, long n_pts,
                               npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_sextic                                          *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a sextic expansion, and then performs diffraction     *
     *      correction on diffracted data using this. Mathematical definitions *
     *      are given in the comment at the start of this file.                *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*(A_2+x2*A_4));
        psi_odd  = kd*x2*x*(A_1+x2*A_3);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

complex double _fresnel_sextic_norm(double* x_arr, char* T_in, double* w_func,
                                    double rcpr_D, double A_0, double A_1,
                                    double A_2, double A_3, double A_4,
                                    double dx, double rcpr_F, double kd,
                                    long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_sextic                                          *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a sextic expansion, and then performs diffraction     *
     *      correction on diffracted data using this. The normalization scheme *
     *      is then applied to the output. Mathematical definitions can be     *
     *      found in the comment at the start of this file.                    *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*(A_2+x2*A_4));
        psi_odd  = kd*x2*x*(A_1+x2*A_3);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += exp_negative_psi+exp_positive_psi;

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double _fresnel_octic(double* x_arr, char* T_in, double* w_func,
                              double rcpr_D, double A_0, double A_1,
                              double A_2, double A_3, double A_4, double A_5,
                              double A_6, double dx, double rcpr_F, double kd,
                              long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_sextic                                          *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a sextic expansion, and then performs diffraction     *
     *      correction on diffracted data using this. Mathematical definitions *
     *      are given in the comment at the start of this file.                *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*(A_2+x2*(A_4+x2*A_6)));
        psi_odd  = kd*x2*x*(A_1+x2*(A_3+x2*A_5));

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

complex double _fresnel_octic_norm(double* x_arr, char* T_in, double* w_func,
                                   double rcpr_D, double A_0, double A_1,
                                   double A_2, double A_3, double A_4,
                                   double A_5, double A_6, double dx,
                                   double rcpr_F, double kd,
                                   long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_sextic                                          *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a sextic expansion, and then performs diffraction     *
     *      correction on diffracted data using this. The normalization scheme *
     *      is then applied to the output. Mathematical definitions can be     *
     *      found in the comment at the start of this file.                    *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        psi_even = kd*x2*(A_0+x2*(A_2+x2*(A_4+x2*A_6)));
        psi_odd  = kd*x2*x*(A_1+x2*(A_3+x2*A_5));

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += exp_negative_psi+exp_positive_psi;

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function, that is where w_func = 1.                            */
    T_out += *(complex double *)T_in;
    norm  += 1.0;

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double _fresnel_transform_newton(double* x_arr, double* phi_arr,
                                         char* T_in, double* w_func, double kD,
                                         double r, double B, double D,
                                         double EPS, long toler, double dx,
                                         double rcpr_F, long n_pts,
                                         npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_newton                                          *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function computes the Fresnel inverse transform using         *
     *      Newton-Raphson to calculate the stationary value of the Fresnel    *
     *      kernel. Mathematical definitions are in the comment at the top     *
     *      of this file.                                                      *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi;
    complex double T_out, exp_psi;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = -(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi     = Fresnel_Psi_Func(kD, r, x_arr[i], phi, phi_arr[i], B, D);
        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_psi * *(complex double *)(T_in + j*T_in_steps);
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

complex double _fresnel_transform_newton_norm(double* x_arr, double* phi_arr,
                                              char* T_in, double* w_func,
                                              double kD, double r, double B,
                                              double D, double EPS, long toler,
                                              double dx, double rcpr_F,
                                              long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Function:                                                              *
     *      _fresnel_transform_newton_norm                                     *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function computes the Fresnel inverse transform using         *
     *      Newton-Raphson to calculate the stationary value of the Fresnel    *
     *      kernel. The normalization scheme is then applied. Mathematical     *
     *      definitions are given in the comment at the start of this file.    *
     **************************************************************************/

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi;
    complex double T_out, exp_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j     = -(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi     = Fresnel_Psi_Func(kD, r, x_arr[i], phi, phi_arr[i], B, D);
        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the norm using a Riemann sum as well.                     */
        norm   += exp_psi;

        /*  Approximate the integral with a Riemann Sum.                      */
        T_out += exp_psi * *(complex double *)(T_in + j*T_in_steps);
        j     += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number.                        */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

/*  End of the Include Guard.                                                 */
#endif