/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H
#define RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H

/*  Various trig functions, complex variables, and more found here.           */
#include <math.h>
#include <complex.h>

/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  This void function takes a pointer as an input. The pointers values       *
 *  are changed within the code, and there is no need to "return" anything.   *
 *  Hence, this function has no "return" statement.                           */
static void get_arr(double* x_arr, double dx, long nw_pts)
{
    /***************************************************************************
     *  Function:                                                              *
     *      get_arr                                                            *
     *  Inputs:                                                                *
     *      x_arr   (Pointer to Double):                                       *
     *          The array that stores the ring radius points, in kilometers.   *
     *      dx      (Double):                                                  *
     *          The spacing between two ring intercept points, in kilometers.  *
     *      nw_pts  (Long):                                                    *
     *          Half of the number of points in a given window, rounded down.  *
     *          For example, for a window with 31 points, nw_pts should be 15. *
     *  Outputs:                                                               *
     *      None. The x_arr pointer is altered within the code.                *
     *  Purpose:                                                               *
     *      This computes an array of length nw_pts, with values ranging from  *
     *      -nw_pts*dx to zero. Due to symmetry in the reconstruction, only    *
     *      the left half of a given window is needed, and hence this returns  *
     *      only half of the array. The values 0 to nw_pts*dx are omitted.     *
     **************************************************************************/
    long i;

    /*  Loop over the input array, and assign the correct ring radius (km)  */
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
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      F           (Double):                                              *
     *          The Fresnel Scale for the center of the window, in kilometers. *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses the classic Fresnel quadratic aproximation to   *
     *      the Fresnel Kernel to perform diffraction correction for the given *
     *      input data. The C99 standard (Or higher) is required, as Complex   *
     *      Variables are used. The transform is defined as:                   *
     *                 W/2                                                     *
     *                  -                                                      *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i pi/2 (r-r_0)^2/F^2) dx_0   *
     *               | |                                                       *
     *                -                                                        *
     *               -W/2                                                      *
     *                                                                         *
     *  Where T_hat is the diffracted data, w is the window function, r is     *
     *  the ring intercept point, and r_0 is a dummy variable of integration.  *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, rcpr_F, rcpr_F2;
    complex double T_out, exp_negative_ix;

    /*  Initialize the T_out variable to zero, so we can loop over later. */
    T_out = 0.0;

    j = -n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprical of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/F;
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.  */
    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computation roughly in half.              */  
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      F           (Double):                                              *
     *          The Fresnel Scale for the center of the window, in kilometers. *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses the classic Fresnel quadratic aproximation to   *
     *      the Fresnel Kernel to perform diffraction correction for the given *
     *      input data. The output is the then normalized by the window width. *
     *      The normalization is given as follows:                             *
     *                                                                         *
     *                |     _ +infinity           |                            *
     *                |    | |                    |                            *
     *                |    |    exp(-i psi(x)) dx |                            *
     *                |  | |                      |                            *
     *                |   -  -infinity            |                            *
     *      Norm =  __________________________________                         *
     *              |    -  +W/2                    |                          *
     *              |   | |                         |                          *
     *              |   |    w(x) exp(-i psi(x)) dx |                          *
     *              | | |                           |                          *
     *              |  -   -W/2                     |                          *
     *                                                                         *
     *      Where psi is the Fresnel Kernel, w is the window function,         *
     *      and W is the window width. |f| denotes the absolute value of f.    *
     *      Variables are used. The transform is defined as:                   *
     *                 W/2                                                     *
     *                  -                                                      *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i pi/2 (r-r_0)^2/F^2) dx_0   *
     *               | |                                                       *
     *                -                                                        *
     *               -W/2                                                      *
     *                                                                         *
     *  Where T_hat is the diffracted data, r is the ring intercept point, and *
     *  r_0 is a dummy variable of integration.                                *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, rcpr_F, rcpr_F2;
    complex double T_out, exp_negative_ix, norm;

    /*  Initialize T_out and norm to zero, so we can loop over later. */
    T_out  = 0.0;
    norm   = 0.0;

    j = -n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprical of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/F;
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        norm  += 2.0*exp_negative_ix;
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;
        norm  += 1.0;

        /*  The integral in the numerator of norm evaluates to F sqrt(2). Use *
         *  this in the calculation of the normalization. The cabs function   *
         *  computes the absolute value of complex number.                    */
        norm = SQRT_2 / cabs(norm);

        /*  Multiply result by the coefficient found in the Fresnel inverse.  *
         *  The 1/F term is omitted, since the F in the norm cancels this.    */
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
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a cubic expansion, and then performs diffraction      *
     *      correction on diffracted data using this. The transform is:        *
     *                 W/2                                                     *
     *                  -                                                      *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x,x_0)) dx_0           *
     *               | |                                                       *
     *                -                                                        *
     *               -W/2                                                      *
     *                                                                         *
     *  Where T_hat is the diffracted data, w is the window function, r is     *
     *  the ring intercept point, and r_0 is a dummy variable of integration.  *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi use Horner's Method for Polynomial Computation.       */
        psi_even = kd*x2*A_0;
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

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *      _fresnel_transform_quartic_norm                                    *
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      A_2         (Double):                                              *
     *          The coefficient of the x^4 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this.The output is the then    *
     *      normalized by the window width. The normalization is given by:     *
     *                                                                         *
     *                |     _ +infinity           |                            *
     *                |    | |                    |                            *
     *                |    |    exp(-i psi(x)) dx |                            *
     *                |  | |                      |                            *
     *                |   -  -infinity            |                            *
     *      Norm =  __________________________________                         *
     *              |    -  +W/2                    |                          *
     *              |   | |                         |                          *
     *              |   |    w(x) exp(-i psi(x)) dx |                          *
     *              | | |                           |                          *
     *              |  -   -W/2                     |                          *
     *                                                                         *
     *      Where psi is the Fresnel Kernel, w is the window function,         *
     *      and W is the window width. |f| denotes the absolute value of f.    *
     *      Variables are used. The transform is defined as:                   *
     *                                                                         *
     *                  -  +W/2                                                *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x)) dx_0               *
     *               | |                                                       *
     *                -   -W/2                                                 *
     *                                                                         *
     *                                                                         *
     *  Where T_hat is the diffracted data, r is the ring intercept point, and *
     *  r_0 is a dummy variable of integration.                                *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
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

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;
        norm  += 1.0;

        /*  The integral in the numerator of norm evaluates to F sqrt(2). Use *
         *  this in the calculation of the normalization. The cabs function   *
         *  computes the absolute value of complex number.                    */
        norm = SQRT_2 / cabs(norm);

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      A_2         (Double):                                              *
     *          The coefficient of the x^4 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this. The transform is:        *
     *                 W/2                                                     *
     *                  -                                                      *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x,x_0)) dx_0           *
     *               | |                                                       *
     *                -                                                        *
     *               -W/2                                                      *
     *                                                                         *
     *  Where T_hat is the diffracted data, w is the window function, r is     *
     *  the ring intercept point, and r_0 is a dummy variable of integration.  *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi use Horner's Method for Polynomial Computation.       */
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

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *      _fresnel_transform_quartic_norm                                    *
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      A_2         (Double):                                              *
     *          The coefficient of the x^4 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this.The output is the then    *
     *      normalized by the window width. The normalization is given by:     *
     *                                                                         *
     *                |     _ +infinity           |                            *
     *                |    | |                    |                            *
     *                |    |    exp(-i psi(x)) dx |                            *
     *                |  | |                      |                            *
     *                |   -  -infinity            |                            *
     *      Norm =  __________________________________                         *
     *              |    -  +W/2                    |                          *
     *              |   | |                         |                          *
     *              |   |    w(x) exp(-i psi(x)) dx |                          *
     *              | | |                           |                          *
     *              |  -   -W/2                     |                          *
     *                                                                         *
     *      Where psi is the Fresnel Kernel, w is the window function,         *
     *      and W is the window width. |f| denotes the absolute value of f.    *
     *      Variables are used. The transform is defined as:                   *
     *                                                                         *
     *                  -  +W/2                                                *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x)) dx_0               *
     *               | |                                                       *
     *                -   -W/2                                                 *
     *                                                                         *
     *                                                                         *
     *  Where T_hat is the diffracted data, r is the ring intercept point, and *
     *  r_0 is a dummy variable of integration.                                *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi use Horner's Method for Polynomial Computation.       */
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

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;
        norm  += 1.0;

        /*  The integral in the numerator of norm evaluates to F sqrt(2). Use *
         *  this in the calculation of the normalization. The cabs function   *
         *  computes the absolute value of complex number.                    */
        norm = SQRT_2 / cabs(norm);

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      A_2         (Double):                                              *
     *          The coefficient of the x^4 term in the expansion for psi.      *
     *      A_3         (Double):                                              *
     *          The coefficient of the x^5 term in the expansion for psi.      *
     *      A_4         (Double):                                              *
     *          The coefficient of the x^6 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this. The transform is:        *
     *                 W/2                                                     *
     *                  -                                                      *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x,x_0)) dx_0           *
     *               | |                                                       *
     *                -                                                        *
     *               -W/2                                                      *
     *                                                                         *
     *  Where T_hat is the diffracted data, w is the window function, r is     *
     *  the ring intercept point, and r_0 is a dummy variable of integration.  *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero so we can loop over later.                   */
    T_out = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi use Horner's Method for Polynomial Computation.       */
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

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
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
     *      _fresnel_transform_quartic_norm                                    *
     *  Inputs:                                                                *
     *      x_arr       (Pointer to Double):                                   *
     *          The ring radii within the given window.                        *
     *      T_in        (Pointer to Char):                                     *
     *          The raw diffraction-limited data that is to be corrected.      *
     *      w_func      (Pointer to Double):                                   *
     *          Pre-computed window function. Should have the same number of   *
     *          points as the x_arr pointer.                                   *
     *      rcpr_D      (Double):                                              *
     *          1/D, where D is the distant from the spacecraft to the ring    *
     *          intercept point, in kilometers.                                *
     *      A_0         (Double):                                              *
     *          The coefficient of the x^2 term in the expansion for psi.      *
     *      A_1         (Double):                                              *
     *          The coefficient of the x^3 term in the expansion for psi.      *
     *      A_2         (Double):                                              *
     *          The coefficient of the x^4 term in the expansion for psi.      *
     *      dx          (Double):                                              *
     *          Spacing between two consecutive ring intercept points, km.     *
     *      rcpr_F      (Double):                                              *
     *          The reciprocal of the Fresnel scale in kilometers.             *
     *      kd          (Double):                                              *
     *          The wavenumber, k, weighted by the spacecraft-ring distance, D.*
     *      n_pts       (Long):                                                *
     *          Half the number of points in the window, rounded down.         *
     *      T_in_steps  (npy_intp (Equivalent to Long)):                       *
     *          The number of steps in memory to get from the nth data point   *
     *          to the (n+1)th data point in the T_in variable (See above).    *
     *  Outputs:                                                               *
     *      T_out       (Complex Double):                                      *
     *          The diffraction corrected profile.                             *
     *  Purpose:                                                               *
     *      This function uses Legendre polynomials to approximate the Fresnel *
     *      kernel up to a quartic expansion, and then performs diffraction    *
     *      correction on diffracted data using this.The output is the then    *
     *      normalized by the window width. The normalization is given by:     *
     *                                                                         *
     *                |     _ +infinity           |                            *
     *                |    | |                    |                            *
     *                |    |    exp(-i psi(x)) dx |                            *
     *                |  | |                      |                            *
     *                |   -  -infinity            |                            *
     *      Norm =  __________________________________                         *
     *              |    -  +W/2                    |                          *
     *              |   | |                         |                          *
     *              |   |    w(x) exp(-i psi(x)) dx |                          *
     *              | | |                           |                          *
     *              |  -   -W/2                     |                          *
     *                                                                         *
     *      Where psi is the Fresnel Kernel, w is the window function,         *
     *      and W is the window width. |f| denotes the absolute value of f.    *
     *      Variables are used. The transform is defined as:                   *
     *                                                                         *
     *                  -  +W/2                                                *
     *                 | |                                                     *
     *      T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(x)) dx_0               *
     *               | |                                                       *
     *                -   -W/2                                                 *
     *                                                                         *
     *                                                                         *
     *  Where T_hat is the diffracted data, r is the ring intercept point, and *
     *  r_0 is a dummy variable of integration.                                *
     **************************************************************************/
    /*  Declare all necessary variables.   */
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    j = -n_pts;

    /*  Use a Riemann Sum to approximation the Fresnel Inverse Integral.      */
    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        /*  Compute psi use Horner's Method for Polynomial Computation.       */
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

        /*  Add the central point in the Riemann sum. This is center of the   *
         *  window function, that is where w_func = 1.                        */
        T_out += *(complex double *)T_in;
        norm  += 1.0;

        /*  The integral in the numerator of norm evaluates to F sqrt(2). Use *
         *  this in the calculation of the normalization. The cabs function   *
         *  computes the absolute value of complex number.                    */
        norm = SQRT_2 / cabs(norm);

        /*  Multiply result by the coefficient found in the Fresnel inverse.  */
        T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

/*  End of the Include Guard. */
#endif