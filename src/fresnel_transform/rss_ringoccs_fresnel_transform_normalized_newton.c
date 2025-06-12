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

/*  Inverse Fresnel transform via Newton-Raphson with window normalization.   */
void
rssringoccs_Fresnel_Transform_Normalized_Newton(
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

    /*  The integrand for the Riemann sum. For the inverse transform this is  *
     *  the product of the tapering function, the transmittance, and the      *
     *  (positive) stationary complex exponentiated Fresnel kernel.           */
    tmpl_ComplexDouble integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. We must compute everything from -W / 2 to +W / 2.          */
    offset = center - ((nw_pts - 1) >> 1);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (n = 0; n < nw_pts; ++n)
    {
        /*  This function does not assume the ideal geometry present in MTR86 *
         *  where u . y = 0, u being the vector from the spacecraft to the    *
         *  ring intercept point. Instead we compute using the full Fresnel   *
         *  kernel. To do this requires rho, rho0, and R as vectors in their  *
         *  Cartesian coordinates. Compute this.                              */
        const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polard(
            tau->rho_km_vals[offset], tau->phi_deg_vals[offset]
        );

        const tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polard(
            tau->rho_km_vals[center], tau->phi_deg_vals[offset]
        );

        const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
            tau->rx_km_vals[offset],
            tau->ry_km_vals[offset],
            tau->rz_km_vals[offset]
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

        /*  The inverse transform is approximated by negating the Fresnel     *
         *  phase, psi. We have:                                              *
         *                                                                    *
         *                   -             _________                          *
         *                  | | ^         / 2 pi     exp(-i (psi_s-pi/4))     *
         *      T(rho) ~=   |   T(r0) r  / --------- -------------------- dr0 *
         *                | |          \/  |psi_s''|     | R - rho_s |        *
         *                 -                                                  *
         *                                                                    *
         *  The kernel for the inverse transform is hence simply the complex  *
         *  conjugate of the kernel for the forward transform. Compute this.  */
        tmpl_CDouble_ConjugateSelf(&ker);

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
/*  End of rssringoccs_Fresnel_Transform_Normalized_Newton.                   */
