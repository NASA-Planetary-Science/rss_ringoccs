/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/
#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

#define RSSRINGOCCS_RCPR_SQRT_TWO (+7.0710678118654752440084436210484903928E-01)

void
rssringoccs_Fresnel_Transform_Legendre_Even_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
)
{
    /*  Declare all necessary variables. n and m are used for indexing.       */
    size_t n;
    size_t m = n_pts;

    /*  k is used for polynomial evaluation using Horner's method. This will  *
     *  index the coefficients.                                               */
    unsigned int k;

    /*  psi is the Fresnel kernel, computed using Legendre polynomials.       */
    double psi;

    /*  This is the multiplicative factor in front of the integral. The total *
     *  factor is (1 + i) / 2F, where F is the Fresnel scale. Since we are    *
     *  normalizing by the window width, the actual scale factor is           *
     *  (1 + i) sqrt(2) / |norm|, where norm is the complex normalization     *
     *  factor obtained by integrating w * exp(-i psi) over a free space      *
     *  region, w being the window function.                                  */
    double scale_factor;

    /*  We will use symmetry to simplify and speed up the calculation. We     *
     *  have psi = L0 x^2 + L1 x^3 + L3 X^4 + ..., replacing x with -x, which *
     *  is equivalent to replacing (r - r0) / D with (r0 - r) / D, leads to   *
     *  psi = L0 x^2 - L1 x^3 + L3 x^4 - ..., the odd terms become negative   *
     *  and the even terms stay the same. We compute the sum of the even      *
     *  terms as psi_even, and the sum of the odd terms as psi_odd. The left  *
     *  side (negative x) is then psi = psi_even - psi_odd, and the right     *
     *  side (positive x) is psi = psi_even + psi_odd. Using such symmetry    *
     *  cuts the size of the for loop in half and saves us some time. These   *
     *  variables are for exp(-i psi(-x)) and exp(-i psi(+x)), respectively.  */
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;

    /*  Dummy variable of integration. This will hold T_hat * w * exp(-i psi).*/
    tmpl_ComplexDouble integrand;

    /*  The complex normalization factor, which is the free space integral.   *
     *  We will compute this using a Riemann sum. Initialize to zero to start.*/
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    const double rcpr_D = 1.0 / tau->D_km_vals[center];

    /*  The factor k * D is used frequently. Save some redundant computations *
     *  by computing this outside of the for loop.                            */
    const double kD = tau->k_vals[center] * tau->D_km_vals[center];

    /*  Initialize T_out to zero so we can loop over it. We will compute      *
     *  T_out, which is equal to the inverse Fresnel transform, via a         *
     *  Riemann sum, summing over all of the data in the current window.      */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Use a Riemann Sum to approximate the inverse Fresnel transform.       */
    for (n = 0; n < n_pts; ++n)
    {
        /*  x_arr stores r - r0 as r varies from r0 - W/2 to r0, and hence is *
         *  always negative. The polynomial is in terms of (r - r0) / D.      *
         *  Compute this, and remember that this value is negative.           */
        const double x = x_arr[n] * rcpr_D;

        /*  We compute the even and odd parts of the expansion separately.    *
         *  Because of this we will use Horner's method in x^2, not x.        */
        const double x2 = x*x;

        /*  The constant and linear terms of the expansion are zero, so there *
         *  is a multiplicative factor of kD * x^2 for psi. Compute this.     */
        const double kD_times_x2 = kD * x2;

        /*  The current diffraction data points. One for the left side of the *
         *  window, and one for the right. These correspond to negative x and *
         *  positive x, respectively.                                         */
        const tmpl_ComplexDouble T_left = tau->T_in[center - m];
        const tmpl_ComplexDouble T_right = tau->T_in[center + m];

        /*  Compute psi using Horner's Method. Initialize the even and odd    *
         *  parts of psi. This function assumes the degree of the expansion   *
         *  is even, meaning the last coefficients is for an even degree term.*
         *  The second to last coefficient is thus for an odd degree term.    *
         *  Initialize accordingly.                                           */
        double psi_even = coeffs[tau->order - 1];
        double psi_odd = coeffs[tau->order - 2];

        /*  Perform Horner's method.                                          */
        for (k = 3; k < tau->order - 1; k += 2)
        {
            psi_even = psi_even*x2 + coeffs[tau->order - k];
            psi_odd = psi_odd*x2  + coeffs[tau->order - k - 1];
        }

        /*  There zeroth coefficient is for the x^2 term which is part of the *
         *  even expansion. Perform the final iteration of Horner's method.   */
        psi_even = psi_even*x2 + coeffs[0];

        /*  Scale by x to finish the computation of psi_odd. That is, we have *
         *                                                                    *
         *                (N-1)/2                                             *
         *                 -----                                              *
         *                 \            2n + 1                                *
         *      psi   =    /      L    x                                      *
         *         odd     -----   2n+1                                       *
         *                 n = 0                                              *
         *                                                                    *
         *                (N-1)/2                                             *
         *                 -----                                              *
         *                 \            2n                                    *
         *            = x  /      L    x                                      *
         *                 -----   2n+1                                       *
         *                 n = 0                                              *
         *                                                                    *
         *  We have thus far computed the bottom sum, but without the factor  *
         *  of x on the outside. Finish the computation, multiply by x.       */
        psi_odd *= x;

        /*  Recalling that x is negative, to compute the left side of the sum *
         *  (where x is negative), we just need to do psi_even - psi_odd. We  *
         *  also need the scale factor kD x^2. Finish the computation for psi.*/
        psi = kD_times_x2 * (psi_even - psi_odd);

        /*  Accounting for the window function, we have w(-x) exp(-i psi(-x)).*
         *  Since window functions are real-valued and non-negative, this is  *
         *  just the polar form of a complex number. Compute from this.       */
        w_exp_minus_psi_left = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  The right side of the window has x positive (r > r0). psi_even is *
         *  not changed by this, but psi_odd flips in sign. Compute.          */
        psi = kD_times_x2 * (psi_even + psi_odd);

        /*  Window functions are symmetric: w(-x) = w(x). We can compute via  *
         *  the polar form method with the same element of the window array.  */
        w_exp_minus_psi_right = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  The complex normalization is the integral over free space. This   *
         *  is computed via a Riemann sum. Add both the left and right sides, *
         *  w(+/- x) exp(-i psi(+/- x)), to the norm factor. The "dx" term is *
         *  factored into the Riemann sum at the end of this for-loop.        */
        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_left);
        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_right);

        /*  The integrand for the left part of the window is:                 *
         *      T_hat(-x) w(-x) exp(-i psi(-x))                               *
         *  The integrand for the right part of the window is similar:        *
         *      T_hat(-x) w(-x) exp(-i psi(-x))                               *
         *  We sum this into T_out to compute the Riemann sum. The "dx"       *
         *  factor is once again ignored for now.                             */
        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_left, T_left);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_right, T_right);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  n is the index for the w_func data, m is the index for the x      *
         *  variable and T_hat. n is incremented, m is decremented since we   *
         *  start at the edge of the window and move towards the center.      */
        m--;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    tmpl_CDouble_AddTo(&tau->T_out[center], &tau->T_in[center]);
    tmpl_CDouble_AddTo_Real(&norm, 1.0);

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    scale_factor = RSSRINGOCCS_RCPR_SQRT_TWO / tmpl_CDouble_Abs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
