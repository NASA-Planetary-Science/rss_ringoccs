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
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform(rssringoccs_TAUObj *tau,
                              const double *x_arr,
                              const double *w_func,
                              size_t n_pts,
                              size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, n;

    /*  Variables for the real and imaginary parts of the output.             */
    double real, imag;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprocal of F as a variable and compute with that.                  */
    const double rcpr_fresnel_scale = 1.0 / tau->F_km_vals[center];

    /*  The forward transform can be computed by negating the psi factor. Set *
     *  the sign to 1 for the inverse transform and -1 for the forward one.   */
    const double sign = (tau->use_fwd ? -1.0 : 1.0);
    const double psi_factor = sign * rcpr_fresnel_scale * rcpr_fresnel_scale;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprocal of F as a variable and compute with that.                  */
    const double scale = 0.5 * tau->dx_km * rcpr_fresnel_scale;

    /*  exp_negative_ix is used for the Fresnel kernel.                       */
    tmpl_ComplexDouble exp_negative_ipsi, integrand;

    /*  Start with the central point in the Riemann sum. This is center of    *
     *  window function. That is, where w_func = 1. This is just T_in at      *
     *  the central point. This also initializes T_out.                       */
    tau->T_out[center] = tau->T_in[center];

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    n = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  The x array passed to us contains (pi/2)(r - r0)^2. The Fresnel   *
         *  approximation is psi = (pi/2) * (r - r0)^2 / F^2.                 */
        const double psi = psi_factor * x_arr[m];

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ipsi = tmpl_CDouble_Polar(w_func[m], -psi);

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        integrand = tmpl_CDouble_Add(tau->T_in[center-n], tau->T_in[center+n]);
        tmpl_CDouble_MultiplyBy(&integrand, &exp_negative_ipsi);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  Decrement the offset index to the next point.                     */
        n--;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    real = scale * (tau->T_out[center].dat[0] - tau->T_out[center].dat[1]);
    imag = scale * (tau->T_out[center].dat[0] + tau->T_out[center].dat[1]);
    tau->T_out[center].dat[0] = real;
    tau->T_out[center].dat[1] = imag;
}
/*  End of rssringoccs_Fresnel_Transform.                                     */
