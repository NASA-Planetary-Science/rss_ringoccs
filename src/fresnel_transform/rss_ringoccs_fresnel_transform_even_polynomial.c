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
 *               rss_ringoccs_fresnel_transform_even_polynomial               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Uses even degree polynomials to approximate the Fresnel kernel at the *
 *      stationary azimuth angle and perform diffraction correction.          *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Fresnel_Transform_Even_Polynomial                         *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using even degree polynomials.        *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          A pointer to a Tau object. This contains all of the geometry and  *
 *          diffraction limited data, and this function will write the newly  *
 *          reconstructed data to the T_out array.                            *
 *      x_arr (const double * TMPL_RESTRICT const):                           *
 *          The array r[n] - r[center], where r is the radius. n = 0 hence    *
 *          corresponds to the left-most edge of the window, n = n_pts + 1    *
 *          represents the center of the window.                              *
 *      w_func (const double * TMPL_RESTRICT const):                          *
 *          The window function, pre-computed across the current window. The  *
 *          value w_func[n] corresponds to the window at x_arr[n] (see above).*
 *      coeffs (const double * TMPL_RESTRICT const):                          *
 *          The coefficients for the polynomial approximation. There must be  *
 *          at least tau->order elements in the array. coeffs[0] represents   *
 *          the constant coefficient, coeffs[tau->order - 1] corresponds to   *
 *          the coefficient of the highest order term.                        *
 *      n_pts (size_t):                                                       *
 *          The number of points in the x_arr and w_func arrays. There are    *
 *          2 * n_pts + 1 points total in the window, n_pts to the left of    *
 *          the center, n_pts to the right, and the center itself.            *
 *      center (size_t):                                                      *
 *          The index for the center of the window. There must be             *
 *          n_pts points to the left and right of the center in the data.     *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_complex.h:                                                       *
 *          tmpl_CDouble_Polar:                                               *
 *              Computes z = r * exp(i theta), with theta in radians.         *
 *          tmpl_CDouble_AddTo:                                               *
 *              Performs z += w, in-place.                                    *
 *          tmpl_CDouble_AddTo_Real:                                          *
 *              Performs z += r, in-place.                                    *
 *          tmpl_CDouble_Multiply:                                            *
 *              Complex multiplication, z = w0 * w1.                          *
 *          tmpl_CDouble_Rect:                                                *
 *              Creates a complex number from two real ones, z = x + iy.      *
 *  Method:                                                                   *
 *      The inverse Fresnel transform is given by:                            *
 *                                                                            *
 *                         infinity                                           *
 *                            -                                               *
 *                   1 + i   | |                                              *
 *          T(rho) = -----   |   T_hat(r0) w(r - r0) exp(-i psi(r, r0)) dr0   *
 *                     2F  | |                                                *
 *                          -                                                 *
 *                        -infinity                                           *
 *                                                                            *
 *      where T_hat is the diffracted data, w is the window function, r is    *
 *      the ring intercept point, and r0 is a dummy variable of integration.  *
 *      psi is the Fresnel Kernel, and exp is the exponential function. The   *
 *      scale factor F is the Fresnel scale which is dependent on the         *
 *      geometry of the occultation.                                          *
 *                                                                            *
 *      This function uses a user-provided polynomial approximation for the   *
 *      Fresnel kernel. It is most commonly used with Legendre polynomials    *
 *      and Chebyshev polynomials of the second kind, but you may use         *
 *      whatever you want.                                                    *
 *                                                                            *
 *      The azimuth angle corresponding to r is the value phi_s such that     *
 *      d psi / d phi = 0, the stationary azimuth angle. The function         *
 *      psi(r, r0) can be approximated as a sum using Legendre polynomials    *
 *      and Chebyshev polynomials of the second kind. We have:                *
 *                                                                            *
 *                                             N                              *
 *                     P (A) - A P   (A)     -----                            *
 *                      n         n+1        \                                *
 *          L (A, B) = ----------------- - B /     P     (A) P (A)            *
 *           n               n + 2           -----  n+k-1     k               *
 *                                           n = 0                            *
 *                                                                            *
 *                                                                            *
 *                     P (A) - A P   (A)      -                     -         *
 *                      n         n+1        |  U   (A) - 2 P   (A)  |        *
 *                   = ----------------- - B |   n+2         n+2     |        *
 *                           n + 2            -                     -         *
 *                                                                            *
 *      where P_n is the nth Legendre polynomial and U_n is the nth           *
 *      Chebyshev polynomial of the second kind. A and B are defined in terms *
 *      of the geometry of the occultation. psi is then approximated by:      *
 *                                                                            *
 *                            N                         n + 2                 *
 *                          -----            -        -                       *
 *                          \               |  r - r0  |                      *
 *          psi(r, r0) = kD /      L (A, B) |  ------  |                      *
 *                          -----   n        -    D   -                       *
 *                          n = 0                                             *
 *                                                                            *
 *      Where D is the distance between the spacecraft and the ring point.    *
 *      We commonly use this function by setting coeffs[n] = L_n(A, B).       *
 *      Note, this is not required. The coeffs array may be arbitrary.        *
 *                                                                            *
 *      The evaluation of the polynomial is done using Horner's method. Since *
 *      (-x)^n = x^n for even n and (-x)^m = -x^m for odd m we may use        *
 *      symmetry to cut the computation roughly in half. That is, we          *
 *      compute the odd part of the polynomial, psi_odd, and the even part,   *
 *      psi_even, as rho varies from rho0 - W/2 to rho0, W being the window   *
 *      width, rho0 being the center of the window. On the left we have:      *
 *                                                                            *
 *          psi = psi_even + psi_odd                                          *
 *                                                                            *
 *      On the right, since the sign of rho - rho0 flips, we have:            *
 *                                                                            *
 *          psi = psi_even - psi_odd                                          *
 *                                                                            *
 *      psi_even and psi_odd are computed using Horner's method for the left  *
 *      side of the window, and then psi is computed using these formulas.    *
 *      The final integral is then computed using a Riemann sum.              *
 *  Notes:                                                                    *
 *      1.) The Legendre approximation assumes the first iteration of         *
 *          Newton's method with guess phi = phi0 is sufficient for finding   *
 *          the stationary azimuth angle, the angle phi_s with                *
 *          d psi / d phi = 0. For low ring opening angles B we require more  *
 *          iterations, often 3 or 4, in order to converge to the root. In    *
 *          these scenarios the Legendre approximation fails.                 *
 *      2.) There are geometries where the quadratic Fresnel approximation is *
 *          better than the higher order Legendre approximations. This is     *
 *          because the quadratic Fresnel approximation is indeed the         *
 *          quadratic term of the true stationary Fresnel kernel, psi with    *
 *          phi = phi_s, but the Legendre approximation only works when the   *
 *          first Newton iterate for phi accurately approximates phi_s.       *
 *      3.) There are no checks for NULL pointers. You are responsible for    *
 *          calling this function with pointers that point to valid data.     *
 *      4.) This function may also be used with the Lagrange interpolating    *
 *          polynomials that are described in the appendix of MTR86.          *
 *      5.) This function does not normalize by the window width, meaning for *
 *          very coarse resolutions the Riemann sum will be roughly zero.     *
 *  References:                                                               *
 *      1.) Maguire, R., French, R. (2025)                                    *
 *          Applications of Legendre Polynomials for Fresnel Inversion        *
 *              and Occultation Observations                                  *
 *                                                                            *
 *          Derivation of the Legendre approximation is given. This function  *
 *          implements some of the mathematics in this paper.                 *
 *                                                                            *
 *      2.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. rss_ringoccs implements many of the       *
 *          ideas found in this article.                                      *
 *                                                                            *
 *      3.) Goodman, J. (2005)                                                *
 *          Introduction to Fourier Optics                                    *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *                                                                            *
 *          Covers most of the theory behind diffraction and the application  *
 *          of Fourier analysis to optics. The Fresnel transform is given an  *
 *          in-depth treatise in this book.                                   *
 *                                                                            *
 *      4.) McQuarrie, Donald (2003),                                         *
 *          Mathematical Methods for Scientists and Engineers                 *
 *          University Science Books, ISBN 1-891389-29-7                      *
 *                                                                            *
 *          Excellent introductory text on mathematical physics. A detailed   *
 *          discussion of orthogonal polynomials can be found in Chapter 14:  *
 *          Orthogonal Polynomials and Sturm-Liouville Problems.              *
 *                                                                            *
 *      5.) Arfken, G., Weber, H., Harris, F. (2013)                          *
 *          Mathematical Methods for Physicists, Seventh Edition              *
 *          Academic Press, Elsevier                                          *
 *                                                                            *
 *          Standard textbook on mathematical physics, often used in          *
 *          introductory graduate courses. See Chapter 18: More Special       *
 *          Functions, Section 4: Chebyshev Polynomials                       *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_config.h:                                                        *
 *          Header file providing TMPL_RESTRICT.                              *
 *  2.) tmpl_complex.h:                                                       *
 *          Complex numbers declared here, as are arithmetic functions.       *
 *  3.) rss_ringoccs_fresnel_transform.h:                                     *
 *          Prototype for the function given here.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2025/04/09 (Ryan Maguire):                                                *
 *      Cleaned up comments, added references, improved organization.         *
 *  2025/04/10 (Ryan Maguire):                                                *
 *      Changed function name to allow for arbitrary polynomials. This still  *
 *      works with the Legendre approximation, but other methods may be used. *
 ******************************************************************************/

/*  TMPL_RESTRICT defined here. This expands to "restrict" if the C99, or     *
 *  higher, standard is supported, and nothing otherwise. This allows         *
 *  rss_ringoccs to be compiled using C89 compilers.                          */
#include <libtmpl/include/tmpl_config.h>

/*  Complex numbers and functions provided here.                              */
#include <libtmpl/include/tmpl_complex.h>

/*  Function prototype provided here.                                         */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Perform diffraction correction using even degree polynomials.             */
void
rssringoccs_Fresnel_Transform_Even_Polynomial(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
)
{
    /*  Declare all necessary variables. n and m are used for indexing.       */
    size_t n;
    size_t m = n_pts;

    /*  k is used for polynomial evaluation using Horner's method. This will  *
     *  index the coefficients.                                               */
    unsigned int k;

    /*  psi is the Fresnel kernel, computed using a polynomial approximation. */
    double psi;

    /*  We will use symmetry to simplify and speed up the calculation. We     *
     *  have psi = L0 x^2 + L1 x^3 + L3 X^4 + ..., replacing x with -x, which *
     *  is equivalent to replacing (r - r0) / D with (r0 - r) / D, leads to   *
     *  psi = L0 x^2 - L1 x^3 + L3 x^4 - ..., the odd terms become negative   *
     *  and the even terms stay the same. We compute the sum of the even      *
     *  terms as psi_even, and the sum of the odd terms as psi_odd. The left  *
     *  side (negative x) is then psi = psi_even + psi_odd, and the right     *
     *  side (positive x) is psi = psi_even - psi_odd. Using such symmetry    *
     *  cuts the size of the for loop in half and saves us some time. These   *
     *  variables are for exp(-i psi(-x)) and exp(-i psi(+x)), respectively.  */
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;

    /*  Dummy variable of integration. This will hold T_hat * w * exp(-i psi).*/
    tmpl_ComplexDouble integrand;

    /*  Multiplicative factor that appears outside of the integral.           */
    const double scale_factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  Division is more expension than division, so store the reciprocal     *
     *  of D as a variable and compute with that.                             */
    const double rcpr_D = 1.0 / tau->D_km_vals[center];

    /*  The factor k * D is used frequently. Save some redundant computations *
     *  by computing this outside of the for loop.                            */
    const double kD = tau->k_vals[center] * tau->D_km_vals[center];

    /*  Initialize T_out to zero so we can loop over it. We will compute      *
     *  T_out, which is equal to the inverse Fresnel transform, via a         *
     *  Riemann sum, summing over all of the data in the current window.      */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Use a Riemann Sum to approximate the inverse Fresnel transform.       */
    for (n = 0; n < n_pts; ++n)
    {
        /*  x_arr stores r - r0 as r varies from r0 - W/2 to r0, and hence is *
         *  always negative. The polynomial is in terms of (r - r0) / D.      *
         *  Compute this, and remember that this value is negative.           */
        const double x = x_arr[n] * rcpr_D;

        /*  We compute the even and odd parts of the expansion separately.    *
         *  Because of this we will use Horner's method in x^2, not x.        */
        const double x2 = x*x;

        /*  The constant and linear terms of the expansion are zero, so there *
         *  is a multiplicative factor of kD * x^2 for psi. Compute this.     */
        const double kD_times_x2 = kD * x2;

        /*  The current diffraction data points. One for the left side of the *
         *  window, and one for the right. These correspond to negative x and *
         *  positive x, respectively.                                         */
        const tmpl_ComplexDouble T_left = tau->T_in[center - m];
        const tmpl_ComplexDouble T_right = tau->T_in[center + m];

        /*  Compute psi using Horner's Method. Initialize the even and odd    *
         *  parts of psi. This function assumes the degree of the expansion   *
         *  is even, meaning the last coefficients is for an even degree term.*
         *  The second to last coefficient is thus for an odd degree term.    *
         *  Initialize accordingly.                                           */
        double psi_even = coeffs[tau->order - 1];
        double psi_odd = coeffs[tau->order - 2];

        /*  Perform Horner's method.                                          */
        for (k = 3; k < tau->order - 1; k += 2)
        {
            psi_even = psi_even*x2 + coeffs[tau->order - k];
            psi_odd = psi_odd*x2 + coeffs[tau->order - k - 1];
        }

        /*  The zeroth coefficient is for the x^2 term which is part of the   *
         *  even expansion. Perform the final iteration of Horner's method.   */
        psi_even = psi_even*x2 + coeffs[0];

        /*  Scale by x to finish the computation of psi_odd. That is, we have *
         *                                                                    *
         *                (N-1)/2                                             *
         *                 -----                                              *
         *                 \            2n + 1                                *
         *      psi   =    /      L    x                                      *
         *         odd     -----   2n+1                                       *
         *                 n = 0                                              *
         *                                                                    *
         *                (N-1)/2                                             *
         *                 -----                                              *
         *                 \            2n                                    *
         *            = x  /      L    x                                      *
         *                 -----   2n+1                                       *
         *                 n = 0                                              *
         *                                                                    *
         *  We have thus far computed the bottom sum, but without the factor  *
         *  of x on the outside. Finish the computation, multiply by x.       */
        psi_odd *= x;

        /*  Recalling that x is negative, to compute the left side of the sum *
         *  (where x is negative), we just need to do psi_even + psi_odd. We  *
         *  also need the scale factor kD x^2. Finish the computation for psi.*/
        psi = kD_times_x2 * (psi_even + psi_odd);

        /*  Accounting for the window function, we have w(-x) exp(-i psi(-x)).*
         *  Since window functions are real-valued and non-negative, this is  *
         *  just the polar form of a complex number. Compute from this.       */
        w_exp_minus_psi_left = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  The right side of the window has x positive (r > r0). psi_even is *
         *  not changed by this, but psi_odd flips in sign. Compute.          */
        psi = kD_times_x2 * (psi_even - psi_odd);

        /*  Window functions are symmetric: w(-x) = w(x). We can compute via  *
         *  the polar form method with the same element of the window array.  */
        w_exp_minus_psi_right = tmpl_CDouble_Polar(w_func[n], -psi);

        /*  The integrand for the left part of the window is:                 *
         *                                                                    *
         *      T_hat(-x) w(-x) exp(-i psi(-x))                               *
         *                                                                    *
         *  The integrand for the right part of the window is similar:        *
         *                                                                    *
         *      T_hat(+x) w(+x) exp(-i psi(+x))                               *
         *                                                                    *
         *  We sum this into T_out to compute the Riemann sum. The "dx"       *
         *  term is contained in the scale_factor variable above.             */
        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_left, T_left);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_right, T_right);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  n is the index for the w_func and x_arr data, m is the index for  *
         *  the T_hat variable. n is incremented, m is decremented since we   *
         *  start at the edge of the window and move towards the center.      */
        m--;
    }

    /*  Add the central point in the Riemann sum. This is the center of the   *
     *  window function. That is, where w_func = 1.                           */
    tmpl_CDouble_AddTo(&tau->T_out[center], &tau->T_in[center]);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
/*  End of rssringoccs_Fresnel_Transform_Even_Polynomial.                     */
