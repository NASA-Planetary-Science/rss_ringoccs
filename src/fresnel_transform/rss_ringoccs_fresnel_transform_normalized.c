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

/*  Complex numbers provided here.                                            */
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform_Normalized(
    rssringoccs_TAUObj * const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
)
{
    /*  Declare all necessary variables. m and n are used for indexing.       */
    size_t m, n;

    /*  The scale factor for the integral, which is 1 / 2F times the          *
     *  magnitude of the complex normalization.                               */
    double scale;

    /*  Variables for the real and imaginary parts of the output.             */
    double real, imag;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprocal of F as a variable and compute with that.                  */
    const double rcpr_fresnel_scale = 1.0 / tau->F_km_vals[center];

    /*  The forward transform can be computed by negating the psi factor. Set *
     *  the sign to 1 for the inverse transform and -1 for the forward one.   */
    const double sign = (tau->use_fwd ? -1.0 : 1.0);
    const double psi_factor = sign * rcpr_fresnel_scale * rcpr_fresnel_scale;

    /*  exp_negative_ipsi is the complex exponentiation Fresnel kernel, and   *
     *  and integrand is a helper variable used for computing the Fresnel     *
     *  transform via a Riemann sum.                                          */
    tmpl_ComplexDouble exp_negative_ipsi, integrand;

    /*  Initialize norm to zero, so we can loop over later.                   */
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;

    /*  At the center of the window we have r = r0, and hence w(r - r0) = 1   *
     *  and psi(r, r0) = 0, meaning exp(i psi) = 1. The contribution to the   *
     *  Riemann sum is just T_in[center]. Initialize the output to this value.*/
    tau->T_out[center] = tau->T_in[center];

    /*  From symmetry we need only compute over the range r0 to r0 + W/2. We  *
     *  start the offset at n_pts and decrement this to zero.                 */
    n = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m < n_pts; ++m)
    {
        /*  The x array passed to us contains (pi/2)(r - r0)^2. The Fresnel   *
         *  approximation is psi = (pi/2) * (r - r0)^2 / F^2.                 */
        const double psi = psi_factor * x_arr[m];

        /*  w(r - r0) exp(-i psi(r, r0)) is in polar form. Compute using this.*/
        exp_negative_ipsi = tmpl_CDouble_Polar(w_func[m], -psi);

        /*  Compute denominator portion of norm using a Riemann Sum. The      *
         *  numerator has a closed form solution of sqrt(2) F. We can skip    *
         *  the explicit computation.                                         */
        tmpl_CDouble_AddTo(&norm, &exp_negative_ipsi);

        /*  Take advantage of the symmetry of the problem. Since (-x)^2 = x^2,*
         *  and since window functions are symmetric, w(-x) = w(x), if we set *
         *  K(r - r0) = w(r - r0) exp(-i (pi / 2) (r - r0)^2 / F^2), we have: *
         *      K(-t) T_hat(-t) + K(t) T_hat(t) = K(t)(T_hat(-t) + T_hat(t))  *
         *  and hence we can compute T_hat(-t) + T_hat(t) first and then      *
         *  scale the result by K(t), where t is the appropriate index for    *
         *  the data points. This cuts the computation roughly in half.       */
        integrand = tmpl_CDouble_Add(tau->T_in[center-n], tau->T_in[center+n]);
        tmpl_CDouble_MultiplyBy(&integrand, &exp_negative_ipsi);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  Decrement the offset index to the next point.                     */
        n--;
    }

    /*  Compute 2 * norm + 1. norm is a complex, do this component-wise. The  *
     *  times 2 factor is because we cut the above window in half, using the  *
     *  symmetry that comes from (-x)^2 = x^2. The plus one comes from the    *
     *  center of the window, which was also skipped in the for-loop above.   */
    norm.dat[0] = 2.0 * norm.dat[0] + 1.0;
    norm.dat[1] *= 2.0;

    /*  Compute the real scale factor, 1 / (sqrt(2) |norm|).                  */
    scale = tmpl_Double_Rcpr_Sqrt_Two / tmpl_CDouble_Abs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1 / F term is omitted, since the F in the norm cancels this.      */
    real = scale * (tau->T_out[center].dat[0] - tau->T_out[center].dat[1]);
    imag = scale * (tau->T_out[center].dat[0] + tau->T_out[center].dat[1]);
    tau->T_out[center].dat[0] = real;
    tau->T_out[center].dat[1] = imag;
}
/*  End of rssringoccs_Fresnel_Transform_Normalized.                          */
