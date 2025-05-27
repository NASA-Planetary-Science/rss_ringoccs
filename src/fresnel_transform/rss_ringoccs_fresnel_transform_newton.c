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
rssringoccs_Fresnel_Transform_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
)
{
    /*  m is used for indexing the Riemann sum, offset is the index for the   *
     *  off-center point in the window, the dummy variable of integration.    */
    size_t m, offset;

    /*  Scale factor for the integral is (1 + i) / 2F, where F is the Fresnel *
     *  scale. The Riemann sum also has a factor of "dx", so the scalar part  *
     *  of the factor is dx / 2F. The "1 + i" part is included later.         */
    const double factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  Complex variable for the express T_hat * w * exp(-i psi). T_hat is    *
     *  the diffracted data, w is the window function, and psi is the         *
     *  Fresnel kernel.                                                       */
    tmpl_ComplexDouble integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - ((n_pts - 1) >> 1);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m < n_pts; ++m)
    {
        /*  This function does not assume the ideal geometry present in MTR86 *
         *  where u . y = 0, u being the vector from the spacecraft to the    *
         *  ring intercept point. Instead we compute using the full Fresnel   *
         *  kernel. To do this requires rho, rho0, and R as vectors in their  *
         *  Cartesian coordinates. Compute this.                              */
        const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polar(
            tau->rho_km_vals[offset],
            tau->phi_deg_vals[offset] * tmpl_Double_Deg_To_Rad
        );

        const tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polar(
            tau->rho_km_vals[center],
            tau->phi_deg_vals[offset] * tmpl_Double_Deg_To_Rad
        );

        const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
            tau->rx_km_vals[offset],
            tau->ry_km_vals[offset],
            tau->rz_km_vals[offset]
        );

        /*  Compute the Fresnel kernel at the stationary azimuth angle, phi_s.*/
        const double psi = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
        );

        /*  The complex exponentiated Fresnel kernel is simply exp(-i psi).   *
         *  This is scaled by the window function and integrated over.        */
        const tmpl_ComplexDouble exp_psi = tmpl_CDouble_Polar(w_func[m], -psi);

        /*  The integrand is the diffracted data times the tapered complex    *
         *  exponentiated Fresnel kernel. Compute and add to the sum.         */
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  We start at the left-most edge of the window and move to the      *
         *  right. We increment "offset" up to center + (n_pts - 1) / 2.      */
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
