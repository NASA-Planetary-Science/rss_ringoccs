#include "_diffraction_correction.h"

/*  All of the kernel functions (psi) found here.                             */
#include "_fresnel_kernel.h"

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
complex double
Fresnel_Transform_Perturbed_Newton_Double(double *x_arr, double *phi_arr,
                                          complex double *T_in,
                                          double *w_func, double kD,
                                          double r, double B, double D,
                                          double EPS, long toler,
                                          double dx, double F, long n_pts,
                                          long center, double perturb[5])
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi, x, poly;
    complex double T_out, exp_psi;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Factor for the polynomial perturbation.                           */
        x = (r-x_arr[i])/D;

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute psi and perturb by the requested polynomial.              */
        psi = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);

        /*  Use Horner's method to compute the polynomial.                    */
        poly = x*perturb[4]+perturb[3];
        poly = poly*x + perturb[2];
        poly = poly*x + perturb[1];
        poly = poly*x + perturb[0];

        poly *= kD;

        psi += poly;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
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

complex double
Fresnel_Transform_Perturbed_Newton_Norm_Double(double *x_arr, double *phi_arr,
                                               complex double *T_in,
                                               double *w_func, double kD,
                                               double r, double B, double D,
                                               double EPS, long toler,
                                               double dx, double F, long n_pts,
                                               long center, double perturb[5])
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, x, poly;
    complex double T_out, exp_psi, norm;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Factor for the polynomial perturbation.                           */
        x = (r-x_arr[i])/D;

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute psi and perturb by the requested polynomial.              */
        psi = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);

        /*  Use Horner's method to compute the polynomial.                    */
        poly = x*perturb[4]+perturb[3];
        poly = poly*x + perturb[2];
        poly = poly*x + perturb[1];
        poly = poly*x + perturb[0];

        poly *= kD;

        psi += poly;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
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
    norm = rssringoccs_Sqrt_Two / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}
