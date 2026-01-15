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

/*  TMPL_RESTRICT macro provided here.                                        */
#include <libtmpl/include/tmpl_config.h>

/*  Double precision complex numbers and routines given here.                 */
#include <libtmpl/include/tmpl_complex.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function prototype / forward declaration found here.                      */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Fresnel transform using Newton's method and a Riemann sum.                */
void
rssringoccs_Fresnel_Transform_Newton_Riemann(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    /*  n is used for indexing the Riemann sum, offset is the index for the   *
     *  off-center point in the window, the dummy variable of integration.    */
    size_t n, offset;

    /*  Scale factor for the window normalization, which is the reciprocal of *
     *  the magnitude of the free space integral across the window.           */
    double scale;

    /*  Variable for the complex Fresnel kernel, part of the integrand.       */
    tmpl_ComplexDouble ker;

    /*  The integrand for the Riemann sum. For the Fresnel transform this is  *
     *  the product of the tapering function, the transmittance, and the      *
     *  stationary complex Fresnel kernel.                                    */
    tmpl_ComplexDouble integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. We must compute everything from -W / 2 to +W / 2.          */
    offset = center - ((nw_pts - 1) >> 1);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (n = 0; n < nw_pts; ++n)
    {
        /*  Compute the full Fresnel kernel, not assuming psi'' is constant.  */
        ker = rssringoccs_Fresnel_Kernel(tau, center, offset);

        /*  Lastly, scale the kernel by the window function.                  */
        tmpl_CDouble_MultiplyBy_Real(&ker, w_func[n]);

        /*  The normalization factor is the free space integral, which is     *
         *  also computed using a Riemann sum.                                */
        if (tau->use_norm)
            tmpl_CDouble_AddTo(&norm, &ker);

        /*  The integrand is the product of the complex data, T_in, and the   *
         *  Fresnel kernel (scaled by the tapering function). Add this to     *
         *  the Riemann sum.                                                  */
        integrand = tmpl_CDouble_Multiply(ker, tau->T_in[offset]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  We are moving left-to-right in the data, increment the offset.    */
        ++offset;
    }

    /*  By Babinet's principle, the free space integral is 1. We are          *
     *  integrating over a finite data set, and have introduced a tapering    *
     *  function. The normalization factor is the magnitude of the free space *
     *  integral across the entire real line divided by the tapered integral  *
     *  across the window, which is hence 1 / | norm |. Compute this.         */
    if (tau->use_norm)
    {
        /*  Scale factor computed using Babinet's principle.                  */
        scale = 1.0 / tmpl_CDouble_Abs(norm);

        /*  Scale the Riemann sum by the normalization factor to conclude.    */
        tmpl_CDouble_MultiplyBy_Real(&tau->T_out[center], scale);
    }

    /*  If the integral is not normalized, we simply need to scale by the     *
     *  "dx" factor that is included in the Riemann sum.                      */
    else
        tmpl_CDouble_MultiplyBy_Real(&tau->T_out[center], tau->dx_km);
}
/*  End of rssringoccs_Fresnel_Transform_Newton_Riemann.                      */
