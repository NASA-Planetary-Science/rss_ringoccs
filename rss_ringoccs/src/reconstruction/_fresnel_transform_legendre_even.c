#include <math.h>
#include <complex.h>
#include <rss_ringoccs_math_constants.h>

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
complex double
Fresnel_Transform_Legendre_Even_Double(double *x_arr, complex double *T_in,
                                       double *w_func, double D,
                                       double *coeffs, double dx, double F,
                                       double kd, long n_pts,
                                       unsigned char order, long center)
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

        /*  Compute psi using Horner's Method for Polynomial Computation.  */
        psi_even = coeffs[order-1];
        psi_odd  = coeffs[order-2];
        for (k=3; k<order-1;){
            psi_even = psi_even*x2 + coeffs[order-k];
            psi_odd  = psi_odd*x2 + coeffs[order-k-1];
            k += 2;
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

complex double
Fresnel_Transform_Legendre_Norm_Even_Double(double *x_arr, complex double *T_in,
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
        psi_odd  = coeffs[order-2];
        for (k=3; k<order-1;){
            psi_even = psi_even*x2 + coeffs[order-k];
            psi_odd  = psi_odd*x2 + coeffs[order-k-1];
            k += 2;
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
