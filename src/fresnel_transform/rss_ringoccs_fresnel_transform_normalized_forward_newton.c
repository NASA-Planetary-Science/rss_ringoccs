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

/*  Double precision complex numbers and routines given here.                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Functions for the Fresnel kernel and Fresnel optics found here.           */
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>

/*  2D and 3D vector typedefs, used for the geometry of the inversion.        */
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function prototype / forward declaration found here.                      */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

#if 0

/*  Forward Fresnel transform via Newton-Raphson with window normalization.   */
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

    /*  Scale factor for the window normalization, which is the reciprocal of *
     *  the magnitude of the free space integral across the window.           */
    double scale;

    /*  The integrand for the Riemann sum. For the forward transform this is  *
     *  the product of the tapering function, the transmittance, and the      *
     *  (positive) stationary complex exponentiated Fresnel kernel.           */
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

        /*  Compute the Fresnel kernel at the stationary azimuth angle, phi_s.*
         *  Note, this function returns the scale factor that is present from *
         *  the stationary phase approximation. The general integral is:      *
         *                                                                    *
         *                          -    -                                    *
         *      ^         sin(B)   | |  | |        exp(i  psi)                *
         *      T(rho0) = ------   |    |   T(rho) ----------- drho           *
         *                lambda | |  | |          | R - rho |                *
         *                        -    -                                      *
         *                                                                    *
         *  Where rho is a vector in this expression. Writing                 *
         *  rho = (r cos(phi), r sin(phi)), drho = r dr dphi, and assuming    *
         *  circular symmetry, T(rho) = T(r), we obtain:                      *
         *                                                                    *
         *                         inf         2 pi                           *
         *                          -           -                             *
         *      ^         sin(B)   | |         | | exp(i  psi)                *
         *      T(rho0) = ------   |  r T(r)   |   ----------- dphi dr        *
         *                lambda | |         | |   | R - rho |                *
         *                        -           -                               *
         *                        0           0                               *
         *                                                                    *
         *  The phi integral is handled via stationary phase, producing:      *
         *                                                                    *
         *                                                                    *
         *        2 pi                                                        *
         *         -                          _________                       *
         *        | | exp(i  psi)            / 2 pi     exp(i (psi_s - pi/4)) *
         *        |   ----------- dphi ~=   / --------- -------------------   *
         *      | |   | R - rho |         \/  |psi_s''|    | R - rho_s |      *
         *       -                                                            *
         *       0                                                            *
         *                                                                    *
         *  Where phi_s is the stationary azimuth angle, psi_s and psi_s''    *
         *  the evaluation of psi and psi'' at phi = phi_s, respectively, and *
         *  where rho_s = (r cos(phi_s), r sin(phi_s)). The stationary        *
         *  cylindrical Fresnel kernel is r times this expression. This       *
         *  function computes this quantity.                                  */
        tmpl_ComplexDouble ker = tmpl_Double_Stationary_Cyl_Fresnel_Kernel(
            tau->k_vals[offset], &rho, &rho0, &R, tau->EPS, tau->toler
        );

        /*  Lastly, scale the kernel by the window function.                  */
        tmpl_CDouble_MultiplyBy_Real(&ker, w_func[n]);

        /*  The normalization factor is the free space integral, which is     *
         *  also computed using a Riemann sum.                                */
        tmpl_CDouble_AddTo(&norm, &ker);

        /*  The integrand is the product of the complex data, T_in, and the   *
         *  Fresnel kernel (scaled by the tapering function). Add this to     *
         *  the Riemann sum.                                                  */
        integrand = tmpl_CDouble_Multiply(ker, tau->T_in[offset]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  We are moving left-to-right in the data, increment the offset.    */
        offset += 1;
    }

    /*  By Babinet's principle, the free space integral is 1. We are          *
     *  integrating over a finite data set, and have introduced a tapering    *
     *  function. The normalization factor is the magnitude of the free space *
     *  integral across the entire real line divided by the tapered integral  *
     *  across the window, which is hence 1 / | norm |. Compute this.         */
    scale = 1.0 / tmpl_CDouble_Abs(norm);

    /*  Scale the Riemann sum by the normalization factor to conclude.        */
    tmpl_CDouble_MultiplyBy_Real(&tau->T_out[center], scale);
}
/*  End of rssringoccs_Fresnel_Transform_Normalized_Forward_Newton.           */

#else

/*  Inverse Fresnel transform via Newton-Raphson with window normalization.   */
void
rssringoccs_Fresnel_Transform_Normalized_Forward_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    size_t n;
    double slope;
    double current_psi, next_psi;
    size_t offset = center - ((nw_pts - 1) >> 1);
    const double factor = 0.5 / tau->F_km_vals[center];
    const double threshold = 1.0E-8;

    tmpl_ComplexDouble integrand, left, right;

    const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polard(
        tau->rho_km_vals[center], tau->phi_deg_vals[center]
    );

    const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
        tau->rx_km_vals[center],
        tau->ry_km_vals[center],
        tau->rz_km_vals[center]
    );

    tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polard(
        tau->rho_km_vals[offset], tau->phi_deg_vals[center]
    );

    current_psi = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
        tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
    );

    tau->T_out[center] = tmpl_CDouble_Zero;

    for (n = 0; n < nw_pts - 1; ++n)
    {
        rho = tmpl_2DDouble_Polard(
            tau->rho_km_vals[offset + 1], tau->phi_deg_vals[center]
        );

        next_psi = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
        );

        slope = (next_psi - current_psi) / tau->dx_km;

        if (tmpl_Double_Abs(slope) < threshold)
        {
            integrand = tmpl_CDouble_Polar(tau->dx_km * w_func[n], current_psi);
            tmpl_CDouble_MultiplyBy(&integrand, &tau->T_in[offset]);
        }

        else
        {
            left = tmpl_CDouble_Expi(current_psi);
            right = tmpl_CDouble_Expi(next_psi);
            integrand = tmpl_CDouble_Subtract(right, left);
            tmpl_CDouble_MultiplyBy_Imag(&integrand, -w_func[n] / slope);
            tmpl_CDouble_MultiplyBy(&integrand, &tau->T_in[offset]);
        }

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        offset += 1;

        current_psi = next_psi;
    }

    integrand = tmpl_CDouble_Rect(factor, -factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
/*  End of rssringoccs_Fresnel_Transform_Normalized_Newton.                   */

#endif
