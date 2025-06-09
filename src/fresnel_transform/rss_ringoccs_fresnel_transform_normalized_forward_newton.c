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
#include <libtmpl/include/tmpl_cyl_fresnel_diffraction.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>

#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform_Normalized_Forward_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    /*  n is used for indexing the Riemann sum, offset is the index for the   *
     *  off-center point in the window, the dummy variable of integration.    */
    size_t n, offset;

    /*  The magnitude of the normalization factor, which is the free space    *
     *  integral, and the scale factor, which is sqrt(2) / |norm|.            */
    double abs_norm, scale;

    /*  The integrand for the Riemann sum. For the forward transform this is  *
     *  the product of the tapering function, the transmittance, and the      *
     *  complex exponentiated Fresnel kernel, exp(i psi).                     */
    tmpl_ComplexDouble integrand;

    /*  This function does not assume the ideal geometry present in MTR86     *
     *  where u . y = 0, u being the vector from the spacecraft to the ring   *
     *  intercept point. Instead we compute using the full Fresnel kernel. To *
     *  do this requires rho, rho0, and R as vectors in their Cartesian       *
     *  coordinates. For the forward transform, R and rho0 are fixed, only    *
     *  rho varies. We may compute these outside of the main for-loop.        */
    const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polard(
        tau->rho_km_vals[center], tau->phi_deg_vals[center]
    );

    const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
        tau->rx_km_vals[center],
        tau->ry_km_vals[center],
        tau->rz_km_vals[center]
    );

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. We must compute everything from -W / 2 to +W / 2.          */
    offset = center - ((nw_pts - 1) >> 1);

    /*  Use a Riemann Sum to approximate the forward Fresnel Integral.        */
    for (n = 0; n < nw_pts; ++n)
    {
        /*  We compute psi at the stationary azimuth angle, the angle phi_s   *
         *  such that d psi / d phi = 0. This is done using Newton's method   *
         *  in the function below. The guess for the root is phi0, so create  *
         *  the vector (r cos(phi0), r sin(phi0)) as the starting point for   *
         *  Newton's method.                                                  */
        const tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polard(
            tau->rho_km_vals[offset], tau->phi_deg_vals[center]
        );

        /*  Compute the Fresnel kernel at the stationary azimuth angle, phi_s.*/
        const double psi = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
        );

        /*  The complex exponentiated Fresnel kernel is simply exp(i psi).    *
         *  This is scaled by the window function and integrated over.        */
        const tmpl_ComplexDouble exp_psi = tmpl_CDouble_Polar(w_func[n], psi);

        /*  Compute the norm using a Riemann sum as well.                     */
        tmpl_CDouble_AddTo(&norm, &exp_psi);

        /*  The integrand is w(r - r0) T_hat[r] exp(i psi(rho, rho)). Compute.*/
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);

        /*  We are performing a simple Riemann sum. Add the integrand.        */
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  We are computing left-to-right, increment to the next point.      */
        offset += 1;
    }

    /*  The integral in the numerator of the scale factor is F sqrt(2),       *
     *  the denominator is |norm| dx, so the ratio is sqrt(2) F / |norm| dx.  *
     *  Since the dx is on the bottom, this will cancel the dx that occurs in *
     *  the Riemann sum for T_out. The complex scale factor for the Fresnel   *
     *  transform (the constant term outside the integral) is (1 - i) / 2 F.  *
     *  We therefore have:                                                    *
     *                                                                        *
     *      1 - i numer      1 - i sqrt(2) F                                  *
     *      ----- ----- dx = ----- --------- dx                               *
     *       2 F  demom       2 F  |norm| dx                                  *
     *                                                                        *
     *                            1 - i                                       *
     *                     = --------------                                   *
     *                       sqrt(2) |norm|                                   *
     *                                                                        *
     *  Compute this and scale the result to finish the calculation.          */
    abs_norm = tmpl_CDouble_Abs(norm);
    scale = tmpl_Double_Rcpr_Sqrt_Two / abs_norm;
    integrand = tmpl_CDouble_Rect(scale, -scale);
    tmpl_CDouble_MultiplyBy(&tau->T_out[center], &integrand);
}
/*  rssringoccs_Fresnel_Transform_Normalized_Forward_Newton.                  */
