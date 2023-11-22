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
rssringoccs_Fresnel_Transform_Newton_dD_dphi(rssringoccs_TAUObj *tau,
                                             const double *w_func,
                                             size_t n_pts,
                                             size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, offset;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, D, factor;
    tmpl_ComplexDouble exp_psi, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;
    factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - (n_pts - 1UL) / 2UL;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_dD_dPhi_Newton(
            tau->k_vals[offset],        /* Wavenumber. */
            tau->rho_km_vals[center],   /* Dummy ring radius. */
            tau->rho_km_vals[offset],   /* Ring radius. */
            tau->phi_deg_vals[offset],  /* Dummy azimuth angle. */
            tau->phi_deg_vals[offset],  /* Ring azimuth angle. */
            tau->B_deg_vals[offset],    /* Ring opening angle. */
            tau->rx_km_vals[center],    /* Cassini x coordinate. */
            tau->ry_km_vals[center],    /* Cassini y coordinate. */
            tau->rz_km_vals[center],    /* Cassini z coordinate. */
            tau->EPS,                   /* Allowed error. */
            tau->toler                  /* Maximum number of iterations. */
        );

        D = tmpl_Double_Cyl_Fresnel_Observer_Distance(
            tau->rho_km_vals[offset],   /* Ring radius. */
            phi,                        /* Stationary azimuth angle. */
            tau->rx_km_vals[center],    /* Cassini x coordinate. */
            tau->ry_km_vals[center],    /* Cassini y coordinate. */
            tau->rz_km_vals[center]     /* Cassini z coordinate. */
        );

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],        /* Wavenumber. */
            tau->rho_km_vals[center],   /* Dummy ring radius. */
            tau->rho_km_vals[offset],   /* Ring radius. */
            phi,                        /* Stationary azimuth angle. */
            tau->phi_deg_vals[offset],  /* Ring azimuth angle. */
            tau->B_deg_vals[center],    /* Ring opening angle. */
            D                           /* Observer distance. */
        );

        exp_psi = tmpl_CDouble_Polar(w_func[m], -psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tau->T_out[center] = tmpl_CDouble_Add(tau->T_out[center], integrand);
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
