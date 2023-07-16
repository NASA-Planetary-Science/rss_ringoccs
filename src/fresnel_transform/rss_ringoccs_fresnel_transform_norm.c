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
rssringoccs_Fresnel_Transform_Norm(rssringoccs_TAUObj *tau,
                                   const double *x_arr,
                                   const double *w_func,
                                   size_t n_pts,
                                   size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, n;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2, abs_norm, real_norm;

    /*  exp_negative_ix is the Fresnel kernel, norm is the normalization.     */
    tmpl_ComplexDouble exp_negative_ix, norm, integrand;

    /*  Start with the central point in the Riemann sum. This is center of    *
     *  window function. That is, where w_func = 1. This is just T_in at      *
     *  the central point. This also initializes T_out.                       */
    tau->T_out[center] = tau->T_in[center];

    /*  Initialize norm to zero, so we can loop over later.                   */
    norm = tmpl_CDouble_Zero;

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    n = n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprical of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/tau->F_km_vals[center];
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        x = x_arr[m]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = tmpl_CDouble_Polar(w_func[m], -x);

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        tmpl_CDouble_AddTo(&norm, &exp_negative_ix);

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        integrand = tmpl_CDouble_Add(tau->T_in[center - n],
                                     tau->T_in[center + n]);
        integrand = tmpl_CDouble_Multiply(exp_negative_ix, integrand);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        n--;
    }

    norm = tmpl_CDouble_Multiply_Real(2.0, norm);
    tmpl_CDouble_AddTo_Real(&norm, 1.0);
    abs_norm = tmpl_CDouble_Abs(norm);
    real_norm = tmpl_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = tmpl_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
