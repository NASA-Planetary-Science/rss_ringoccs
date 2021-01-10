/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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

#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
Fresnel_Transform_Newton_Norm_Double(rssringoccs_TAUObj *tau,
                                     double *w_func,
                                     unsigned long n_pts,
                                     unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long m, offset;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, cos_psi, sin_psi, real_norm, abs_norm;
    rssringoccs_ComplexDouble exp_psi, norm, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = rssringoccs_CDouble_Zero;
    norm = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - (long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset],
            tau->phi_rad_vals[offset],
            tau->phi_rad_vals[offset],
            tau->B_rad_vals[center],
            tau->D_km_vals[center],
            tau->EPS,
            tau->toler
        );

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = rssringoccs_Double_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset],
            phi,
            tau->phi_rad_vals[offset],
            tau->B_rad_vals[center],
            tau->D_km_vals[center]
        );

        cos_psi = w_func[m]*rssringoccs_Double_Cos(psi);
        sin_psi = w_func[m]*rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);

        /*  Compute the norm using a Riemann sum as well.                     */
        norm = rssringoccs_CDouble_Add(norm, exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tau->T_out[center] = rssringoccs_CDouble_Add(tau->T_out[center],
                                                     integrand);
        offset += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = rssringoccs_CDouble_Abs(norm);
    real_norm = rssringoccs_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = rssringoccs_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    tau->T_out[center] = rssringoccs_CDouble_Multiply(integrand,
                                                      tau->T_out[center]);
}
