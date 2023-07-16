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
#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform_Legendre_Even_Norm(rssringoccs_TAUObj *tau,
                                                 const double *x_arr,
                                                 const double *w_func,
                                                 const double *coeffs,
                                                 size_t n_pts,
                                                 size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t i, j;
    unsigned int k;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double x, x2, psi;
    double psi_even, psi_odd;
    double abs_norm, real_norm;
    tmpl_ComplexDouble exp_negative_psi, exp_positive_psi;
    tmpl_ComplexDouble norm, integrand;

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    const double rcpr_D = 1.0/tau->D_km_vals[center];

    /*  Initialize T_out to zero so we can loop over later.                   */
    tau->T_out[center] = tmpl_CDouble_Zero;
    norm = tmpl_CDouble_Zero;
    j = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i)
    {
        x  = x_arr[i]*rcpr_D;
        x2 = x*x;

        /*  Compute psi using Horner's Method for Polynomial Computation.  */
        psi_even = coeffs[tau->order-1];
        psi_odd  = coeffs[tau->order-2];
        for (k = 3U; k < tau->order-1; k += 2U)
        {
            psi_even = psi_even*x2 + coeffs[tau->order-k];
            psi_odd  = psi_odd*x2  + coeffs[tau->order-k-1];
        }

        /*  The leading term is x^2, so multiply by this and kD.              */
        psi_even  = psi_even*x2 + coeffs[0];
        psi_even *= tau->k_vals[center]*tau->D_km_vals[center] * x2;
        psi_odd  *= tau->k_vals[center]*tau->D_km_vals[center] * x2 * x;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = tmpl_CDouble_Polar(w_func[i], -psi);

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = tmpl_CDouble_Polar(w_func[i], -psi);

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        tmpl_CDouble_AddTo(&norm, &exp_negative_psi);
        tmpl_CDouble_AddTo(&norm, &exp_positive_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = tmpl_CDouble_Multiply(exp_negative_psi, tau->T_in[center-j]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(exp_positive_psi, tau->T_in[center+j]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        j--;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    tmpl_CDouble_AddTo(&tau->T_out[center], &tau->T_in[center]);
    tmpl_CDouble_AddTo_Real(&norm, 1.0);

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = tmpl_CDouble_Abs(norm);
    real_norm = tmpl_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
