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
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

rssringoccs_ComplexDouble
Fresnel_Transform_Norm_Double(double *x_arr, rssringoccs_ComplexDouble *T_in,
                              double *w_func, double F,
                              unsigned long n_pts, unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long m, n;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2, cos_x, sin_x, abs_norm, real_norm;

    /*  exp_negative_ix is the Fresnel kernel, norm is the normalization.     */
    rssringoccs_ComplexDouble T_out, exp_negative_ix, norm, integrand, arg;

    /*  Initialize T_out and norm to zero, so we can loop over later.         */
    T_out = rssringoccs_CDouble_Zero;
    norm  = rssringoccs_CDouble_Zero;

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    n = n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprical of F as a variable and compute with that.                  */
    rcpr_F  = 1.0/F;
    rcpr_F2 = rcpr_F*rcpr_F;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        x = x_arr[m]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        cos_x = rssringoccs_Double_Cos(x);
        sin_x = rssringoccs_Double_Sin(x);
        arg = rssringoccs_CDouble_Rect(cos_x, -sin_x);
        exp_negative_ix = rssringoccs_CDouble_Multiply_Real(w_func[m], arg);

        /*  Compute denominator portion of norm using a Riemann Sum.          */
        arg  = rssringoccs_CDouble_Multiply_Real(2.0, exp_negative_ix);
        norm = rssringoccs_CDouble_Add(norm, arg);

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        integrand = rssringoccs_CDouble_Add(T_in[center - n], T_in[center + n]);
        integrand = rssringoccs_CDouble_Multiply(exp_negative_ix, integrand);
        T_out = rssringoccs_CDouble_Add(T_out, integrand);
        n -= 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    T_out = rssringoccs_CDouble_Add(T_out, T_in[center]);
    norm  = rssringoccs_CDouble_Add_Real(1.0, norm);
    abs_norm = rssringoccs_CDouble_Abs(norm);
    real_norm = rssringoccs_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    arg   = rssringoccs_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    T_out = rssringoccs_CDouble_Multiply(arg, T_out);
    return T_out;
}
