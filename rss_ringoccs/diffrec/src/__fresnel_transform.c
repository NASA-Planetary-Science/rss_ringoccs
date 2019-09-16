#include "__diffraction_functions.h"

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