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
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Computes the inverse Fresnel transform using a quadratic              *
 *      approximation to the stationary Fresnel kernel.                       *
 *                                                                            *
 *      This is computed by calculating the stationary Fresnel kernel         *
 *      at the extreme end points of the window (rho - w/2 and rho + w/2,     *
 *      where w is the window width and rho is the radius of the point)       *
 *      using Newton-Raphson's method.                                        *
 *                                                                            *
 *      The stationary Fresnel kernel evaluates to zero at the center of the  *
 *      window, meaning we have three points and can compute the unique       *
 *      quadratic polynomial as a function of rho that interpolates these.    *
 *                                                                            *
 *      This is an order of magnitude faster than computing the stationary    *
 *      Fresnel kernel at every point of the window, but is also not very     *
 *      accurate for large windows. For example, Rev007 1km resolution gives  *
 *      a near perfect reconstruction, whereas Rev133 1km resolution is poor. *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a Tau object containing the geometry, power, and     *
 *          phase information needed for reconstruction.                      *
 *      w_func (const double *):                                              *
 *          An array representing the window function for the data set.       *
 *      n_pts (size_t):                                                       *
 *          The number of points in the window.                               *
 *      center (size_t):                                                      *
 *          The index corresponding to the center of the window.              *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       May 4, 2021                                                   *
 ******************************************************************************/

/*  Complex routines provided by libtmpl.                                     */
#include <libtmpl/include/tmpl.h>

/*  Function prototype given here.                                            */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Performs Fresnel inversion by approximation the kernel with a quadratic.  */
void
rssringoccs_Fresnel_Transform_Newton_Quadratic(rssringoccs_TAUObj *tau,
                                               const double *w_func,
                                               size_t n_pts,
                                               size_t center)
{
    /*  Declare all necessary variables. n is used for indexing.              */
    size_t n;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    const size_t start = center - (n_pts - 1UL) / 2UL;
    const size_t end = start + n_pts - (size_t)1;
    size_t offset = start;

    /*  Interpolating coefficients and the values of the kernel.              */
    double C[2], psi_n[2];

    /*  Dummy variable for integrating over and the spacecraft distance.      */
    double x, D;

    /*  Variable for the interpolated Fresnel kernel, and azimuth angle.      */
    double psi, phi;

    /*  Variables for the complex Fresnel integral and integrand being summed.*/
    tmpl_ComplexDouble exp_psi, integrand;

    /*  Scale factors for the interpolation are given by the window size.     */
    const double rcpr_w = 1.0 / tau->w_km_vals[center];
    const double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Scale factor for the Fresnel integral, (1 + i) dx / 2F. The complex   *
     *  part will be computed later, save dx / 2F as a real variable.         */
    const double factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  For quadratic interpolation we use the two extreme endpoints. The     *
     *  center value of psi is exactly zero, giving us three data points.     */
    const double rho[2] = {
        tau->rho_km_vals[center] - 0.5*tau->w_km_vals[center],
        tau->rho_km_vals[center] + 0.5*tau->w_km_vals[center]
    };

    /*  Linear interpolation to compute the corresponding azimuth angles.     *
     *  Compute the slope of rho vs phi.                                      */
    const double num = tau->phi_deg_vals[end] - tau->phi_deg_vals[start];
    const double den = tau->rho_km_vals[end] - tau->rho_km_vals[start];
    const double slope = num / den;

    /*  Compute phi using the slope-intercept formula of the line.            */
    const double phi0[4] = {
        (rho[0] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset],
        (rho[1] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset]
    };

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Compute the values of the Fresnel kernel at the two extreme values.   */
    for (n = (size_t)0; n < (size_t)2; ++n)
    {
        /*  Compute the stationary azimuth angle using Newton-Raphson.        */
        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_D_Newton(
            tau->k_vals[center],        /*  Wavenumber.                       */
            tau->rho_km_vals[center],   /*  Ring radius.                      */
            rho[n],                     /*  Dummy radius.                     */
            phi0[n],                    /*  Initial stationary azimuth guess. */
            phi0[n],                    /*  Dummy azimuthal angle.            */
            tau->B_deg_vals[center],    /*  Ring opening angle.               */
            tau->rx_km_vals[center],    /*  x-coordinate of spacecraft.       */
            tau->ry_km_vals[center],    /*  y-coordinate of spacecraft.       */
            tau->rz_km_vals[center],    /*  z-coordinate of spacecraft.       */
            tau->EPS,                   /*  Epsilon error for Newton's method.*/
            tau->toler                  /*  Maximum number of iterations.     */
        );

        /*  The spacecraft distance is computable from the x, y, z values.    */
        D = tmpl_Double_Cyl_Fresnel_Observer_Distance(
            rho[n],                     /* Ring radius.                       */
            phi,                        /* Stationary azimuth angle.          */
            tau->rx_km_vals[center],    /* Cassini x coordinate.              */
            tau->ry_km_vals[center],    /* Cassini y coordinate.              */
            tau->rz_km_vals[center]     /* Cassini z coordinate.              */
        );

        /*  The Fresnel kernel with variables in radians and kilometers.      */
        psi_n[n] = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],        /*  Wavenumber.                       */
            tau->rho_km_vals[center],   /*  Ring radius.                      */
            rho[n],                     /*  Dummy radius.                     */
            phi,                        /*  Stationary azimuthal angle.       */
            phi0[n],                    /*  Dummy azimuthal angle.            */
            tau->B_deg_vals[center],    /*  Ring opening angle.               */
            D                           /*  Ring-Spacecraft distance.         */
        );
    }

    /*  Compute the coefficients for the interpolating quadratic.             */
    C[0] = (psi_n[1] - psi_n[0]) * 0.5 * rcpr_w;
    C[1] = (psi_n[1] + psi_n[0]) * 0.5 * rcpr_w_sq;

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.       */
    for (n = (size_t)0; n < n_pts; ++n)
    {
        /*  Variable being integrated is rho - rho0, compute this.            */
        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];

        /*  Horner's method to evaluate the quadratic interpolation.          */
        psi = x*(C[0] + x*C[1]);

        /*  Compute w(rho - rho0) * e^{-i psi}, the Fresnel kernel scaled     *
         *  by the tapering (window) function.                                */
        exp_psi = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  The integrand is the transmittance T_hat times the kernel.        */
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);

        /*  Riemann sum, add this current data point to the output.           */
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  And lastly move on to the next data point.                        */
        offset++;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
/*  End of rssringoccs_Fresnel_Transform_Quadratic.                           */
