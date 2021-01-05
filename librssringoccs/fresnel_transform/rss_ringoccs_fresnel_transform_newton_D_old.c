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
Fresnel_Transform_Newton_D_Old_Double(double *x_arr, double *phi_arr,
                                      rssringoccs_ComplexDouble *T_in,
                                      double *w_func, double k, double r,
                                      double B, double EPS, unsigned long toler,
                                      double dx, double F, unsigned long n_pts,
                                      unsigned long center, double rx,
                                      double ry, double rz)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, j;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, x, y, drx, dry, D, exp_psi_re, exp_psi_im;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;
    double rcpr_F = 1.0/F;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i)
    {
        x = x_arr[i] * cos(phi);
        y = x_arr[i] * sin(phi);
        dx = x-rx;
        dy = y-ry;
        D = sqrt(dx*dx + dy*dy + rz*rz);
        kD = k*D;

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi_D_Old(k, r, x_arr[i], phi_arr[i],
                                               phi_arr[i], B, EPS, toler,
                                               rx, ry, rz);

        x = x_arr[i] * cos(phi);
        y = x_arr[i] * sin(phi);
        dx = x-rx;
        dy = y-ry;
        D = sqrt(dx*dx + dy*dy + rz*rz);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = kD*rssringoccs_Double_Fresnel_Psi(1.0, r, x_arr[i], phi,
                                                phi_arr[i], B, D);
        exp_psi_re = cos(psi)*w_func[i];
        exp_psi_im = -sin(psi)*w_func[i];
        exp_psi = rssringoccs_CDouble_Rect(exp_psi_re, exp_psi_im);


        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[j]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        j += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(0.5*dx*rcpr_F, 0.5*dx*rcpr_F);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
