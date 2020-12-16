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
Fresnel_Transform_Perturbed_Newton_Norm_Double(double *x_arr, double *phi_arr,
                                               rssringoccs_ComplexDouble *T_in,
                                               double *w_func, double kD,
                                               double r, double B, double D,
                                               double EPS, unsigned long toler,
                                               unsigned long n_pts,
                                               unsigned long center,
                                               double perturb[5])
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long m, n;

    /*  The Fresnel kernel and the stationary ring azimuth angle.             */
    double psi, phi, x, poly, cos_psi, sin_psi, abs_norm, real_norm;
    rssringoccs_ComplexDouble T_out, exp_psi, norm, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    norm  = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    n = center-(long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  Factor for the polynomial perturbation.                           */
        x = (r-x_arr[m])/D;

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[m], phi_arr[m],
                                         phi_arr[m], B, D, EPS, toler);

        /*  Compute psi and perturb by the requested polynomial.              */
        psi = rssringoccs_Double_Fresnel_Psi(kD, r, x_arr[m], phi,
                                             phi_arr[m], B, D);

        /*  Use Horner's method to compute the polynomial.                    */
        poly  = x*perturb[4]+perturb[3];
        poly  = poly*x + perturb[2];
        poly  = poly*x + perturb[1];
        poly  = poly*x + perturb[0];
        poly *= kD;
        psi  += poly;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[m], exp_psi);

        /*  Compute the norm using a Riemann sum as well.                     */
        norm = rssringoccs_CDouble_Add(norm, exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[m]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        n += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm  = rssringoccs_CDouble_Add_Real(1.0, norm);
    abs_norm = rssringoccs_CDouble_Abs(norm);
    real_norm = rssringoccs_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = rssringoccs_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
