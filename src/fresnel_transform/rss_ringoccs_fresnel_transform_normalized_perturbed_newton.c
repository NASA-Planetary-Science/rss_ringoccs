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
rssringoccs_Fresnel_Transform_Normalized_Perturbed_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, offset;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, x, poly, abs_norm, real_norm;
    tmpl_ComplexDouble exp_psi, norm, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;
    norm = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - (n_pts - 1UL) / 2UL;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m < n_pts; ++m)
    {
        /*  Factor for the polynomial perturbation.                           */
        x = (tau->rho_km_vals[center]-tau->rho_km_vals[offset]) /
            tau->D_km_vals[center];

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = tmpl_Double_Ideal_Stationary_Cyl_Fresnel_Phi_Newton_Deg(
            tau->k_vals[center],        /* Wavenumber. */
            tau->rho_km_vals[center],   /* Dummy radius. */
            tau->rho_km_vals[offset],   /* Ring radius. */
            tau->phi_deg_vals[offset],  /* Dummy azimuthal angle. */
            tau->phi_deg_vals[offset],  /* Ring azimuth angle. */
            tau->B_deg_vals[center],    /* Ring opening angle. */
            tau->D_km_vals[center],     /* Observer distance. */
            tau->EPS,                   /* Allowed error. */
            tau->toler                  /* Max number of iterations. */
        );

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = tmpl_Double_Ideal_Cyl_Fresnel_Psi_Deg(
            tau->k_vals[center],        /* Wavenumber. */
            tau->rho_km_vals[center],   /* Dummy radius. */
            tau->rho_km_vals[offset],   /* Ring radius. */
            phi,                        /* Stationary azimuth angle. */
            tau->phi_deg_vals[offset],  /* Ring azimuth angle. */
            tau->B_deg_vals[center],    /* Ring opening angle. */
            tau->D_km_vals[center]      /* Observer distance. */
        );

        /*  Use Horner's method to compute the polynomial.                    */
        poly  = x*tau->perturb[4] + tau->perturb[3];
        poly  = poly*x + tau->perturb[2];
        poly  = poly*x + tau->perturb[1];
        poly  = poly*x + tau->perturb[0];
        poly *= tau->k_vals[center] * tau->D_km_vals[center];
        psi  += poly;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        exp_psi = tmpl_CDouble_Polar(w_func[m], -psi);

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
