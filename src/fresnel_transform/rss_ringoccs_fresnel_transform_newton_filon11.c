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
 *             rss_ringoccs_fresnel_transform_newton_linear_filon             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Computes the Fresnel transform using a linear Filon-like quadrature.  *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Fresnel_Transform_Newton_Filon11                          *
 *  Purpose:                                                                  *
 *      Performs the Fresnel transform using the Newton-Raphson method to     *
 *      compute the stationary Fresnel phase, and a linear Filon-like         *
 *      quadrature method for the integral.                                   *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const tau):                   *
 *          The Tau object with all of the geometry and diffraction data.     *
 *      w_func (const double * TMPL_RESTRICT const):                          *
 *          The array for the window / tapering function.                     *
 *      center (const size_t):                                                *
 *          The index for the center of the window currently being considered.*
 *          For the forward transform this is the index for "rho0", and for   *
 *          the inverse transform this is the index for "rho".                *
 *      offset (const size_t):                                                *
 *          The index for the variable of integration which runs across the   *
 *          window. This is "rho" for the forward transform and "rho0" for    *
 *          the inverse transform.                                            *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *  Method:                                                                   *
 *      The Fresnel transform is defined as follows:                          *
 *                                                                            *
 *                                                   -                -       *
 *                     -                            |     -      pi -  |      *
 *                    | |                       exp |  i | psi - --- | |      *
 *                    |              ----------     |     -   s   4 -  |      *
 *          ^         |             /   2 pi         -                -       *
 *          T(r ) =   |   r T(r)   / ---------- ------------------------ dr   *
 *             0    | |          \/  | psi '' |       || R - rho  ||          *
 *                   -                    s                    s              *
 *                                                                            *
 *      where psi is the Fresnel phase:                                       *
 *                                                                            *
 *                   -                                          -             *
 *                  |                    R - rho0                |            *
 *          psi = k | || R - rho || - -------------- . (R - rho) |            *
 *                  |                 || R - rho0 ||             |            *
 *                   -                                          -             *
 *                                                                            *
 *      where k is the wavenumber, R is the position vector for the observer, *
 *      rho0 is the where the line of sight from the observer intersects the  *
 *      ring plane, and rho is the dummy variable of integration. The r and   *
 *      r0 variables are given by || rho || and || rho0 ||, respectively.     *
 *      Lastly, rho_s is the stationary point for psi, the point satisfying:  *
 *                                                                            *
 *          d psi                                                             *
 *          ----- = 0                                                         *
 *          d phi                                                             *
 *                                                                            *
 *      where phi is the (dummy) azimuth angle. psi_s is the value of psi at  *
 *      rho = rho_s.                                                          *
 *                                                                            *
 *      The inverse transform is approximated by swapping the rho and rho0    *
 *      variables (hence swapping r and r0 as well) and then integrated the   *
 *      complex conjugate of the integrand:                                   *
 *                                                                            *
 *                                                    -                -      *
 *                     -                             |     -      pi -  |     *
 *                    | |                        exp | -i | psi - --- | |     *
 *                    |               ----------     |     -   s   4 -  |     *
 *                    |     ^        /   2 pi         -                -      *
 *          T(r) ~=   |   r T(r )   / ---------- ------------------------ dr  *
 *                  | |        0  \/  | psi '' |       || R - rho  ||       0 *
 *                   -                     s                     s            *
 *                                                                            *
 *      psi_s is computed using Newton-Raphson with initial guess rho = rho0. *
 *                                                                            *
 *      For large windows (large W), the highly oscillatory nature of the     *
 *      complex exponentiated Fresnel phase means we are unable to use the    *
 *      usual methods of numerical integration, such as the trapezoidal rule  *
 *      or Simpson's rule. Instead, the integral is done using a modified     *
 *      Filon quadrature by doing a linear interpolation of psi and T across  *
 *      each bin. That is, across an interval [Ln, Rn], we interpolate the    *
 *      weight in the integrand (which is everything except the exponential)  *
 *      by a r + b, and interpolate psi by c r + d. We may then approximate   *
 *      the Fresnel transform as follows:                                     *
 *                                                                            *
 *                      R                                                     *
 *                       n                                                    *
 *                      -                                                     *
 *                     | |                 -             -                    *
 *          ^          |   -       -      |    -       -  |                   *
 *          T (r ) =   |  | a r + b | exp | i | c r + d | | dr                *
 *           n  0    | |   -       -      |    -       -  |                   *
 *                    -                    -             -                    *
 *                    L                                                       *
 *                     n                                                      *
 *                                                                            *
 *                     N                                                      *
 *                   -----                                                    *
 *          ^        \     ^                                                  *
 *          T(r ) ~= /     T (r )                                             *
 *             0     -----  n  0                                              *
 *                   n = 0                                                    *
 *                                                                            *
 *      Where N is the number of bins in the window. This mimics the standard *
 *      Filon quadrature method, but we are allowing the slope of the Fresnel *
 *      phase ("c" in the previous equation) to vary from bin to bin. This    *
 *      integral can be evaluated exactly. We have:                           *
 *                                                                            *
 *                      R                                                     *
 *                       n                                                    *
 *                      -                                                     *
 *                     | |                 -             -                    *
 *          ^          |   -       -      |    -       -  |                   *
 *          T (r ) =   |  | a r + b | exp | i | c r + d | | dr                *
 *           n  0    | |   -       -      |    -       -  |                   *
 *                    -                    -             -                    *
 *                    L                                                       *
 *                     n                                                      *
 *                         -                                        -         *
 *                        |      R                     R             |        *
 *                        |       n                     n            |        *
 *                        |      -                     -             |        *
 *                        |     | |                   | |            |        *
 *                    i d |     |     i c r           |   i c r      |        *
 *                 = e    | a   |  r e      d r + b   |  e      d r  |        *
 *                        |   | |                   | |              |        *
 *                        |    -                     -               |        *
 *                        |    L                     L               |        *
 *                        |     n                     n              |        *
 *                         -                                        -         *
 *                                                                            *
 *                           -                                        -       *
 *                      i d |   i R c                i L c             |      *
 *                 = a e    |  e      (1 - i R c) - e     (1 - i L c)  |      *
 *                           -                                        -       *
 *                                      -                  -                  *
 *                             b   i d |    i R c    i L c  |                 *
 *                          + --- e    |  e       - e       |                 *
 *                            i c       -                  _                  *
 *                                                                            *
 *      By definition, we have R c + d = psi_R and L c + d = psi_L. The       *
 *      above simplifies and becomes:                                         *
 *                                                                            *
 *                        -      -                     -      -               *
 *          ^            | i psi  |                   | i psi  |              *
 *          T (r ) = exp |      R | (1 - i R c) - exp |      L | (1 - i L c)  *
 *           n            -      -                     -      -               *
 *                                                                            *
 *                             -                                 -            *
 *                            |       -      -          -      -  |           *
 *                         b  |      | i psi  |        | i psi  | |           *
 *                      + --- |  exp |      R |  - exp |      L | |           *
 *                        i c |       -      -          -      -  |           *
 *                             -                                 -            *
 *                                                                            *
 *      Next, we note that 1 / i = -i, and furthermore that, by definition,   *
 *      we have a R + b = T_R and a L + b = T_L. The previous expression      *
 *      simplified even further and becomes:                                  *
 *                                                                            *
 *                             -          -              -         -          *
 *          ^         i psi   |  a         |    i psi  |  a         |         *
 *          T (r ) =       R  |  - - i T   | -       L |  - - i T   |         *
 *           n  0    e        |  c      R  |   e       |  c      L  |         *
 *                             -          -             -          -          *
 *                   ------------------------------------------------         *
 *                                           c                                *
 *                                                                            *
 *      This ratio, a / c, is the ratio of the slopes for psi and T. The dx   *
 *      term cancels and we have:                                             *
 *                                                                            *
 *              T    -   T                                                    *
 *          a    R        L                                                   *
 *          - = -----------                                                   *
 *          c   psi  - psi                                                    *
 *                 R      L                                                   *
 *                                                                            *
 *      Summing the previous expression over all n would approximate the      *
 *      Fresnel transform provided we have infinite precision. On real        *
 *      computers, using the previous expression when "c" is small results in *
 *      catastrophic cancellation, and a poor numerical integral. For bins    *
 *      where psi_R - psi_L is small (less than 1 / 4 radians), we use the    *
 *      trapezoidal rule to evaluate the integral.                            *
 *  Notes:                                                                    *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. The Fresnel kernel is described here.     *
 *                                                                            *
 *      2.) Goodman, J. (2005)                                                *
 *          Introduction to Fourier Optics                                    *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *                                                                            *
 *          Covers most of the theory behind diffraction and the application  *
 *          of Fourier analysis to optics. The Fresnel transform is given an  *
 *          in-depth treatise in this book.                                   *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       November 12, 2025                                             *
 ******************************************************************************/

/*  TMPL_RESTRICT macro provided here.                                        */
#include <libtmpl/include/tmpl_config.h>

/*  Double precision complex numbers and routines given here.                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Numerical integration tools found here.                                   */
#include <libtmpl/include/tmpl_integration.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>

/*  Function prototype / forward declaration found here.                      */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Inverse Fresnel transform via Newton-Raphson with window normalization.   */
void
rssringoccs_Fresnel_Transform_Newton_Filon11(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    const size_t nw_pts,
    const size_t center
)
{
    double weight, scale;
    double left_psi, right_psi;

    tmpl_ComplexDouble left, right, integrand;

    size_t offset = center - (nw_pts >> 1);
    size_t n;

    tau->T_out[center] = tmpl_CDouble_Zero;

    rssringoccs_Fresnel_Phase_And_Weight(
        tau, center, offset, &weight, &left_psi
    );

    scale = weight * w_func[0];
    left = tmpl_CDouble_Multiply_Real(scale, tau->T_in[offset]);

    for (n = 0; n < nw_pts - 1; ++n)
    {
        rssringoccs_Fresnel_Phase_And_Weight(
            tau, center, offset + 1, &weight, &right_psi
        );

        scale = weight * w_func[n + 1];
        right = tmpl_CDouble_Multiply_Real(scale, tau->T_in[offset + 1]);

        integrand = tmpl_CDouble_Filon11_Integrand(
            left, right, left_psi, right_psi, tau->dx_km
        );

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        left_psi = right_psi;
        left = right;
        ++offset;
    }
}
/*  End of rssringoccs_Fresnel_Transform_Newton_Filon11.                      */
