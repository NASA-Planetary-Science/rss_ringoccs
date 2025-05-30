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
 *      Fresnel_Legendre_Double / Fresnel_Legendre_Double_Norm                *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using Legendre polynomials to         *
 *      approximate the Fresnel kernel. Do to the nature of Legendre          *
 *      polynomials and the first iteration of the Newton-Raphson scheme      *
 *      applied to the Fresnel kernel, this is often extremely accurate and   *
 *      extremely fast.                                                       *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed and rho is the dummy variable of integration  *
 *          which varies from rho0-W/2 to rho0, W being the window width.     *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr. This should *
 *          contain n_pts (see below) number of elements.                     *
 *      D (double):                                                           *
 *          The distance from the spacecraft to the ring-intercept point.     *
 *      coeffs (double *):                                                    *
 *          The coefficients of the polynomial approximation of psi. This     *
 *          should conter order (see below) number of elements.               *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      F (double):                                                           *
 *          The Fresnel scale, in the same units as D and dx.                 *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      order (int):                                                          *
 *          Defined as degree-1, where degree is the degree of the polynomial *
 *          approximation for psi. Order is the highest Legendre polynomial   *
 *          to be used. This is also the size of coeffs (See above).          *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
void
rssringoccs_Fresnel_Transform_Legendre_Odd(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
)
{
    /*  Declare all necessary variables. n, m, and k are used for indexing.   */
    size_t n;
    size_t m = n_pts;
    unsigned int k;

    /*  Variables for the Fresnel kernel and ring radii.                      */
    double psi;
    tmpl_ComplexDouble exp_negative_psi, exp_positive_psi, integrand;

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    const double rcpr_D = 1.0 / tau->D_km_vals[center];
    const double factor = 0.5 * tau->dx_km / tau->F_km_vals[center];
    const double kD = tau->k_vals[center]*tau->D_km_vals[center];

    /*  Initialize T_out to zero so we can loop over later.                   */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (n = 0; n < n_pts; ++n)
    {
        const double x = x_arr[n]*rcpr_D;
        const double x2 = x*x;
        const double kD_times_x2 = kD * x2;

        /*  Compute psi using Horner's Method for Polynomial Computation.     */
        double psi_odd = coeffs[tau->order - 1];
        double psi_even = coeffs[tau->order - 2];

        for (k = 3; k < tau->order - 1; k += 2)
        {
            psi_odd = psi_odd*x2  + coeffs[tau->order - k];
            psi_even = psi_even*x2 + coeffs[tau->order - k - 1];
        }

        /*  The leading term is x^2, so multiply by this and kD.              */
        psi_even *= kD_times_x2;
        psi_odd *= kD_times_x2 * x;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = psi_even - psi_odd;
        exp_negative_psi = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  Compute the right side of exp(-ipsi) using Euler's Formula.       */
        psi = psi_even + psi_odd;
        exp_positive_psi = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = tmpl_CDouble_Multiply(exp_negative_psi, tau->T_in[center-m]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(exp_positive_psi, tau->T_in[center+m]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        m--;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    tmpl_CDouble_AddTo(&tau->T_out[center], &tau->T_in[center]);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
