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
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform_Newton_Norm(rssringoccs_TAUObj *tau,
                                          const double *w_func,
                                          size_t n_pts,
                                          size_t center)
{
    /*  m is used for indexing the Riemann sum, offset is the index for the   *
     *  off-center point in the window, the dummy variable of integration.    */
    size_t m, offset;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double real_norm, abs_norm;
    tmpl_ComplexDouble integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;
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

        /*  Compute the norm using a Riemann sum as well.                     */
        tmpl_CDouble_AddTo(&norm, &exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        offset += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = tmpl_CDouble_Abs(norm);
    real_norm = tmpl_Double_Rcpr_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = tmpl_CDouble_Rect(real_norm, real_norm);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
