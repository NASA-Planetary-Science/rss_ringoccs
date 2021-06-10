#include <math.h>
#include <libtmpl/include/tmpl_complex.h>
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

void
rssringoccs_Fresnel_Transform_Quadratic(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        unsigned long n_pts,
                                        unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, ind[4], offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[2], factor, rcpr_w, rcpr_w_sq, psi_n[4], x;
    double psi, phi, sin_psi, cos_psi;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    tmpl_ComplexDouble exp_psi, integrand;

    rcpr_w = 1.0 / tau->w_km_vals[center];
    rcpr_w_sq = rcpr_w * rcpr_w;
    factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    ind[0] = 0;
    ind[1] = (n_pts-1)/4;
    ind[2] = 3*(n_pts-1)/4;
    ind[3] = n_pts - 1;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - (n_pts - 1UL) / 2UL;

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i < 4; ++i)
    {

        phi = rssringoccs_Newton_Raphson_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            tau->D_km_vals[center],
            tau->EPS,
            tau->toler
        );

        psi_n[i] = rssringoccs_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            phi,
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            tau->D_km_vals[center]
        );
    }

    psi_half_mean = (psi_n[1] + psi_n[2]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[3]) / 2;
    psi_half_diff = psi_n[1] - psi_n[2];
    psi_full_diff = psi_n[0] - psi_n[3];

    C[0] = (8*psi_half_diff - psi_full_diff) * rcpr_w * 0.333333333333333333333;
    C[1] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;

    for (i = 0; i<n_pts; ++i)
    {
        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];
        psi = C[1]*x + C[0];
        psi = psi*x;

        cos_psi = w_func[i]*cos(psi);
        sin_psi = w_func[i]*sin(psi);
        exp_psi = tmpl_CDouble_Rect(cos_psi, -sin_psi);
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tau->T_out[center] = tmpl_CDouble_Add(tau->T_out[center], integrand);
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
