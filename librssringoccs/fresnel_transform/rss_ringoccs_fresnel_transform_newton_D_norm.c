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

rssringoccs_ComplexDouble
Fresnel_Transform_Newton_D_Norm_Double(
    double *x_arr, double *phi_arr, rssringoccs_ComplexDouble *T_in,
    double *w_func, double k, double r, double B, double EPS, long toler,
    long n_pts, long center, double rx, double ry, double rz
){

    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, x, y, dx, dy, D, exp_psi_re, exp_psi_im, abs_norm, real_norm;
    rssringoccs_ComplexDouble T_out, exp_psi, norm, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    norm  = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi_D(k, r, x_arr[i], phi_arr[i],
                                           phi_arr[i], B, EPS, toler,
                                           rx, ry, rz);

        x = x_arr[i] * cos(phi);
        y = x_arr[i] * sin(phi);
        dx = x-rx;
        dy = y-ry;
        D = sqrt(dx*dx + dy*dy + rz*rz);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = rssringoccs_Double_Fresnel_Psi(k, r, x_arr[i], phi,
                                             phi_arr[i], B, D);
        exp_psi_re = cos(psi)*w_func[i];
        exp_psi_im = -sin(psi)*w_func[i];
        exp_psi = rssringoccs_CDouble_Rect(exp_psi_re, exp_psi_im);

        /*  Compute the norm using a Riemann sum as well.                     */
        norm = rssringoccs_CDouble_Add(norm, exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[j]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        j += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = rssringoccs_CDouble_Abs(norm);
    real_norm = rssringoccs_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = rssringoccs_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
