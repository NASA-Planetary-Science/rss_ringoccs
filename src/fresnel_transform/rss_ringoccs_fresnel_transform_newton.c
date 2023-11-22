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

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double*
 *  Purpose:                                                                  *
 *      Performs diffraction correction by using the Newton-Raphson method of *
 *      root-finding to compute the stationary value of the Fresnel kernel.   *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed and rho is the dummy variable of integration  *
 *          which varies from rho0-W/2 to rho+W/2, W being the window width.  *
 *          This should contain n_pts (see below) number of elements.         *
 *      phi_arr (double *):                                                   *
 *          The ring azimuth angle corresponding to rho, in radians.          *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *          This should contain n_pts (see below) number of elements.         *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      r (double):                                                           *
 *          The radius of the current point be reconstructed.                 *
 *      B (double):                                                           *
 *          The ring opening angle, in radians.                               *
 *      D (double):                                                           *
 *          The distance from the spacecraft to the ring-intercept point.     *
 *      EPS (double):                                                         *
 *          The allowable error in the computation of the stationary value of *
 *          the Fresnel kernel.                                               *
 *      toler (long):                                                         *
 *          The maximum number of iterations the Newton-Raphson method is     *
 *          allowed to undergo.                                               *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      F (double):                                                           *
 *          The Fresnel scale, in the same units as D and dx.                 *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
void
rssringoccs_Fresnel_Transform_Newton(rssringoccs_TAUObj *tau,
                                     const double *w_func,
                                     size_t n_pts,
                                     size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi, factor;
    tmpl_ComplexDouble exp_psi, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;
    factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center - (n_pts - 1UL)/2UL;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_Newton(
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
        psi = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],        /* Wavenumber. */
            tau->rho_km_vals[center],   /* Dummy radius. */
            tau->rho_km_vals[offset],   /* Ring radius. */
            phi,                        /* Stationary azimuth angle. */
            tau->phi_deg_vals[offset],  /* Ring azimuth angle. */
            tau->B_deg_vals[center],    /* Ring opening angle. */
            tau->D_km_vals[center]      /* Observer distance. */
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
