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

/*  TMPL_RESTRICT macro provided here.                                        */
#include <libtmpl/include/tmpl_config.h>

/*  Abs function provided here.                                               */
#include <libtmpl/include/tmpl_math.h>

/*  Double precision complex numbers and routines given here.                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Functions for the Fresnel kernel and Fresnel optics found here.           */
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function prototype / forward declaration found here.                      */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Inverse Fresnel transform via Newton-Raphson with window normalization.   */
void
rssringoccs_Fresnel_Transform_Newton_Linear_Filon(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    const double rcpr_dx = 1.0 / tau->dx_km;
    const double threshold = 0.25 * rcpr_dx;
    size_t n, offset;

    double weight, scale;
    double left_psi, right_psi, slope;

    tmpl_ComplexDouble left, right, left_exp_ipsi, right_exp_ipsi;
    tmpl_ComplexDouble left_in, right_in, integrand, complex_slope;

    tau->T_out[center] = tmpl_CDouble_Zero;

    offset = center - ((nw_pts - 1) >> 1);

    rssringoccs_Fresnel_Phase_And_Weight(
        tau, offset, center, &weight, &left_psi
    );

    scale = weight * w_func[0];
    left_exp_ipsi = tmpl_CDouble_Expi(left_psi);

    left_in = tau->T_in[offset];

    tmpl_CDouble_MultiplyBy_Real(&left_in, scale);

    for (n = 0; n < nw_pts - 1; ++n)
    {
        rssringoccs_Fresnel_Phase_And_Weight(
            tau, offset + 1, center, &weight, &right_psi
        );

        scale = weight * w_func[n+1];
        right_exp_ipsi = tmpl_CDouble_Expi(right_psi);

        right_in = tau->T_in[offset + 1];

        tmpl_CDouble_MultiplyBy_Real(&right_in, scale);

        slope = (right_psi - left_psi) * rcpr_dx;

        if (tmpl_Double_Abs(slope) < threshold)
        {
            left = tmpl_CDouble_Multiply(left_in, left_exp_ipsi);
            right = tmpl_CDouble_Multiply(right_in, right_exp_ipsi);
            integrand = tmpl_CDouble_Add(left, right);
            tmpl_CDouble_MultiplyBy_Real(&integrand, 0.5 * tau->dx_km);
        }

        else
        {
            const double rcpr_slope = 1.0 / slope;

            complex_slope = tmpl_CDouble_Subtract(right_in, left_in);
            tmpl_CDouble_MultiplyBy_Real(&complex_slope, rcpr_dx * rcpr_slope);

            left = tmpl_CDouble_Multiply_Imag(1.0, left_in);
            left = tmpl_CDouble_Subtract(left, complex_slope);
            tmpl_CDouble_MultiplyBy_Real(&left, rcpr_slope);
            tmpl_CDouble_MultiplyBy(&left, &left_exp_ipsi);

            right = tmpl_CDouble_Multiply_Imag(1.0, right_in);
            right = tmpl_CDouble_Subtract(right, complex_slope);
            tmpl_CDouble_MultiplyBy_Real(&right, rcpr_slope);
            tmpl_CDouble_MultiplyBy(&right, &right_exp_ipsi);

            integrand = tmpl_CDouble_Subtract(left, right);
        }

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        left_psi = right_psi;
        left_exp_ipsi = right_exp_ipsi;
        left_in = right_in;
        ++offset;
    }
}
/*  End of rssringoccs_Fresnel_Transform_Newton_Linear_Filon.                 */
