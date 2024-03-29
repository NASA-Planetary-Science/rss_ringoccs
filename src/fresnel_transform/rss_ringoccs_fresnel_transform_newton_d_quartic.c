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
rssringoccs_Fresnel_Transform_Newton_D_Quartic(rssringoccs_TAUObj *tau,
                                               const double *w_func,
                                               size_t n_pts,
                                               size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t i, ind[4], offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[4], factor, rcpr_w, rcpr_w_sq, psi_n[4], D, x;
    double psi, phi;
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

        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_dD_dPhi_Newton(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            tau->rx_km_vals[center],
            tau->ry_km_vals[center],
            tau->rz_km_vals[center],
            tau->EPS,
            tau->toler
        );

        D = tmpl_Double_Cyl_Fresnel_Observer_Distance(
            tau->rx_km_vals[center],    /* Ring radius. */
            phi,                        /* Stationary azimuth angle. */
            tau->rx_km_vals[center],    /* Cassini x coordinate. */
            tau->ry_km_vals[center],    /* Cassini y coordinate. */
            tau->rz_km_vals[center]     /* Cassini z coordinate. */
        );

        psi_n[i] = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            phi,
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            D
        );
    }

    psi_half_mean = (psi_n[1] + psi_n[2]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[3]) / 2;
    psi_half_diff = psi_n[1] - psi_n[2];
    psi_full_diff = psi_n[0] - psi_n[3];

    C[0] = (8*psi_half_diff - psi_full_diff) * rcpr_w * 0.333333333333333333333;
    C[1] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[2] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;
    C[3] = (psi_full_mean-4.0*psi_half_mean)*rcpr_w_sq*rcpr_w_sq*21.33333333333;

    for (i = 0; i<n_pts; ++i){

        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];
        psi = C[3];
        psi = psi*x + C[2];
        psi = psi*x + C[1];
        psi = psi*x + C[0];
        psi = psi*x;

        exp_psi = tmpl_CDouble_Polar(w_func[i], -psi);
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tau->T_out[center] = tmpl_CDouble_Add(tau->T_out[center], integrand);
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
