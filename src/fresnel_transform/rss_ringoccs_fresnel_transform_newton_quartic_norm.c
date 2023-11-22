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

#define ONE_THIRD (0.3333333333333333333333333333)
#define FOUR_THIRDS (1.33333333333333333333333333)
#define SIXTEEN_THIRDS (5.3333333333333333333333333)
#define SIXTY_FOUR_THIRDS (21.3333333333333333333333333)

void
rssringoccs_Fresnel_Transform_Newton_Quartic_Norm(rssringoccs_TAUObj *tau,
                                                  const double *w_func,
                                                  size_t n_pts,
                                                  size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t i, ind[4], offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[4], rho[4], abs_norm, real_norm;
    double psi_n[4], psi_half_diff, psi_full_diff;
    double psi, phi;
    double psi_half_mean, psi_full_mean, x;
    tmpl_ComplexDouble exp_psi, norm, integrand;

    const double rcpr_w = 1.0 / tau->w_km_vals[center];
    const double rcpr_w_sq = rcpr_w * rcpr_w;
    const double rcpr_w_cb = rcpr_w_sq * rcpr_w;
    const double rcpr_w_qr = rcpr_w_sq * rcpr_w_sq;

    const size_t quarter = (n_pts - 1) >> 2U;
    ind[0] = 0U;
    ind[3] = n_pts - 1U;

    /*  Special case when the input is one more than a multiple of 4.         */
    if ((n_pts - 1U) & 3U)
    {
        ind[1] = 1U + quarter;
        ind[2] = 1U + 3U*quarter;
    }
    else
    {
        ind[1] = quarter;
        ind[2] = 3U*quarter;
    }

    rho[0] = tau->rho_km_vals[center] - 0.5*tau->w_km_vals[center];
    rho[1] = tau->rho_km_vals[center] - 0.25*tau->w_km_vals[center];
    rho[2] = tau->rho_km_vals[center] + 0.25*tau->w_km_vals[center];
    rho[3] = tau->rho_km_vals[center] + 0.5*tau->w_km_vals[center];

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;
    norm = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - ((n_pts - 1U) >> 1U);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.       */
    for (i = 0; i < 4; ++i)
    {
        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_Newton(
            tau->k_vals[center],                /* Wavenumber. */
            tau->rho_km_vals[center],           /* Dummy radius. */
            rho[i],                             /* Ring radius. */
            tau->phi_deg_vals[offset + ind[i]], /* Dummy azimuthal angle. */
            tau->phi_deg_vals[offset + ind[i]], /* Ring azimuthal angle. */
            tau->B_deg_vals[center],            /* Ring opening angle. */
            tau->D_km_vals[center],             /* Observer distance. */
            tau->EPS,                           /* Allowed error. */
            tau->toler                          /* Max number of iterations. */
        );

        psi_n[i] = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],                /* Wavenumber. */
            tau->rho_km_vals[center],           /* Dummy radius. */
            rho[i],                             /* Ring radius. */
            phi,                                /* Stationary azimuth. */
            tau->phi_deg_vals[offset + ind[i]], /* Ring azimuth. */
            tau->B_deg_vals[center],            /* Ring opening. */
            tau->D_km_vals[center]              /* Observer distance. */
        );
    }

    psi_half_mean = (psi_n[1] + psi_n[2]) * 0.5;
    psi_full_mean = (psi_n[0] + psi_n[3]) * 0.5;
    psi_half_diff = psi_n[2] - psi_n[1];
    psi_full_diff = psi_n[3] - psi_n[0];

    C[0] = (8.0*psi_half_diff - psi_full_diff) * rcpr_w * ONE_THIRD;
    C[1] = (16.0*psi_half_mean - psi_full_mean) * rcpr_w_sq * FOUR_THIRDS;
    C[2] = (psi_full_diff - 2.0*psi_half_diff) * rcpr_w_cb * SIXTEEN_THIRDS;
    C[3] = (psi_full_mean - 4.0*psi_half_mean) * rcpr_w_qr * SIXTY_FOUR_THIRDS;

    for (i = 0; i<n_pts; ++i)
    {
        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];
        psi = x*(C[0] + x*(C[1] + x*(C[2] + x*C[3])));

        exp_psi = tmpl_CDouble_Polar(w_func[i], -psi);
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        tmpl_CDouble_AddTo(&norm, &exp_psi);
        offset += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = tmpl_CDouble_Abs(norm);
    real_norm = tmpl_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = tmpl_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}

#undef ONE_THIRD
#undef FOUR_THIRDS
#undef SIXTEEN_THIRDS
#undef SIXTY_FOUR_THIRDS
