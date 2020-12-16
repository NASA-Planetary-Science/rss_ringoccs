#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stdlib.h>

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Quadratic/Cubic/Quartic                             *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the polynomial interpolation      *
 *      methods described in MTR86. The polynomials are computed by first     *
 *      performing the Newton-Raphson method of root finding for psi, and     *
 *      fitting a quartic polynomial to the end and mid-points of this        *
 *      function across the given window.                                     *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          The ring radius as it varies from rho +/- W/2, W being the window *
 *          width and rho being the center of the window.                     *
 *      phi_arr (double *):                                                   *
 *          The ring azimuth angle corresponding to rho, in radians.          *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. An error check is performed at the *
 *          level above this function, but not explicitly here.               *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr-rho.         *
 *          This should contain n_pts (see below) number of elements.         *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      r (double):                                                           *
 *          The radius of the current point be reconstructed (i.e rho).       *
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
 *          The number of points in the window width. Roughly equal to W/dx.  *
 *      center (long):                                                        *
 *          The index of the center of the diffracted data, T_in. Diffraction *
 *          correction is computed by performing an integral from             *
 *          center-n_pts/2 to center+n_pts/2. If T_in does not have data at   *
 *          these indices, a segmentation fault occurs.                       *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile at the point center (see above).*
 *  Notes:                                                                    *
 *      There are two options for each polynomial, with and without norm. The *
 *      norm options scales the output by a normalization scheme that is a    *
 *      function of the window width. For poor resolutions (10 km or higher)  *
 *      the integral is performed over a very small window, and thus the      *
 *      integral evaluates to roughly zero, and not one. To account for this  *
 *      a factor that is dependent only on the window size (and not the data) *
 *      may be introduced to appropriately scale this back to one. This is    *
 *      default option selected in the Python code.                           *
 ******************************************************************************/

rssringoccs_ComplexDouble
Fresnel_Transform_Quadratic_Double(double *x_arr, double *phi_arr,
                                   rssringoccs_ComplexDouble *T_in,
                                   double *w_func, double kD, double r,
                                   double B, double D,
                                   double EPS, unsigned long toler, double dx,
                                   double F, unsigned long n_pts,
                                   unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C_2, rcpr_F, rcpr_w, rcpr_w_sq;
    double* psi_n;
    double psi, phi, width, sin_psi, cos_psi;
    double psi_half_mean, psi_full_mean;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;

    rcpr_F = 1.0/F;

    /*  Allocate memory for psi_n. Per C99 recommendations, we do not cast    *
     *  malloc since void is safely promoted to whatever type psi_n is.       */
    psi_n = malloc(sizeof(*psi_n) * n_pts);

    width = x_arr[n_pts-1] - x_arr[0];
    rcpr_w = 1.0 / width;
    rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = rssringoccs_Double_Fresnel_Psi(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;

    free(psi_n);

    C_2 = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;

    for (i = 0; i<n_pts; ++i){

        psi = (r-x_arr[i])*(r-x_arr[i])*C_2;
        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[i], exp_psi);
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[j]);
        T_out = rssringoccs_CDouble_Add(T_out, integrand);
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(0.5*dx*rcpr_F, 0.5*dx*rcpr_F);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}

