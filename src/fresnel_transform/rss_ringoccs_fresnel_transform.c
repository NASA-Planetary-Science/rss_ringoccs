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
 *                       rss_ringoccs_fresnel_transform                       *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Uses the classic quadratic approximation to compute the inverse       *
 *      Fresnel transform using a Riemann sum.                                *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Fresnel_Transform                                         *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using the Fresnel approximation.      *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          A pointer to a Tau object. This contains all of the geometry and  *
 *          diffraction limited data, and this function will write the newly  *
 *          reconstructed data to the T_out array.                            *
 *      x_arr (const double * TMPL_RESTRICT const):                           *
 *          The array (pi / 2)(r[n] - r[center])^2, where r is the radius.    *
 *          n = 0 hence corresponds to the left-most edge of the window, and  *
 *          n = n_pts represents the center of the window. The length of the  *
 *          array is given by n_pts, so x_arr holds only the left half of the *
 *          window.                                                           *
 *      w_func (const double * TMPL_RESTRICT const):                          *
 *          The window function, pre-computed across the current window. The  *
 *          value w_func[n] corresponds to the window at x_arr[n] (see above).*
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
 *          tmpl_CDouble_Add:                                                 *
 *              Performs z = w0 + w1.                                         *
 *          tmpl_CDouble_AddTo:                                               *
 *              Performs z += w, in-place.                                    *
 *          tmpl_CDouble_MultiplyBy:                                          *
 *              Performs z *= w, in-place.                                    *
 *  Method:                                                                   *
 *      The inverse Fresnel transform is given by:                            *
 *                                                                            *
 *                           infinity                                         *
 *                              -                                             *
 *                     1 + i   | |                                            *
 *          T_out(r) = -----   |   T_in(r0) w(r - r0) exp(-i psi(r, r0)) dr0  *
 *                       2F  | |                                              *
 *                            -                                               *
 *                          -infinity                                         *
 *                                                                            *
 *      where T_in is the diffracted data, w is the window function, r is     *
 *      the ring intercept point, and r0 is a dummy variable of integration.  *
 *      psi is the Fresnel Kernel, and exp is the exponential function. The   *
 *      scale factor F is the Fresnel scale which is dependent on the         *
 *      geometry of the occultation.                                          *
 *                                                                            *
 *      In ideal scenarios, the Fresnel kernel may be approximated by a       *
 *      simple quadratic:                                                     *
 *                                                                            *
 *                     -      - 2                                             *
 *                pi  | r - r0 |                                              *
 *          psi = --- | ------ |                                              *
 *                 2  |    F   |                                              *
 *                     -      -                                               *
 *                                                                            *
 *      The above integral is then computed via a Riemann sum using this new  *
 *      expression for psi.                                                   *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          occultation accurately in less than a second at 1km resolution.   *
 *      3.) There are no checks for NULL pointers. You are responsible for    *
 *          calling this function with pointers that point to valid data.     *
 *      4.) This function does not normalize by the window width, meaning for *
 *          very coarse resolutions the Riemann sum will be roughly zero.     *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. rss_ringoccs implements many of the       *
 *          ideas found in this article.                                      *
 *                                                                            *
 *      2.) Goodman, J. (2005)                                                *
 *          Introduction to Fourier Optics                                    *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *                                                                            *
 *          Covers most of the theory behind diffraction and the application  *
 *          of Fourier analysis to optics. The Fresnel transform is given an  *
 *          in-depth treatise in this book.                                   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_config.h:                                                        *
 *          Header file providing TMPL_RESTRICT.                              *
 *  2.) tmpl_complex.h:                                                       *
 *          Complex numbers declared here, as are arithmetic functions.       *
 *  3.) rss_ringoccs_tau.h:                                                   *
 *          Tau object typedef provided here.                                 *
 *  4.) rss_ringoccs_fresnel_transform.h:                                     *
 *          Prototype for the function given here.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       May 4, 2021                                                   *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2025/05/29 (Ryan Maguire):                                                *
 *      Cleaned up, added comments / docstring, fixed typos.                  *
 ******************************************************************************/

/*  TMPL_RESTRICT defined here. This expands to "restrict" if the C99, or     *
 *  higher, standard is supported, and nothing otherwise. This allows         *
 *  rss_ringoccs to be compiled using C89 compilers.                          */
#include <libtmpl/include/tmpl_config.h>

/*  Complex numbers and functions provided here.                              */
#include <libtmpl/include/tmpl_complex.h>

/*  Typedef for rssringoccs_TAUObj found here.                               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function prototype provided here.                                         */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Perform the quadratic Fresnel transform on a data set.                    */
void
rssringoccs_Fresnel_Transform(rssringoccs_TAUObj * TMPL_RESTRICT const tau,
                              const double * TMPL_RESTRICT const x_arr,
                              const double * TMPL_RESTRICT const w_func,
                              size_t n_pts,
                              size_t center)
{
    /*  Declare all necessary variables. m and n are used for indexing.       */
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

    /*  The scale factor for the integral is (1 + i) / 2F, and the Riemann    *
     *  sum contains a dx factor, so the final scale is (1 + i) dx / 2F.      */
    const double scale = 0.5 * tau->dx_km * rcpr_fresnel_scale;

    /*  exp_negative_ix is used for the Fresnel kernel.                       */
    tmpl_ComplexDouble exp_negative_ipsi, integrand;

    /*  Start with the central point in the Riemann sum. This is the center   *
     *  of the window function, that is, where w_func = 1. This is just T_in  *
     *  at the central point, which we can use to initialize T_out.           */
    tau->T_out[center] = tau->T_in[center];

    /*  From symmetry we need only compute from W/2 to zero. Start at n_pts,  *
     *  the edge of the window, and work your way in towards the center.      */
    n = n_pts;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m < n_pts; ++m)
    {
        /*  The x array passed to us contains (pi/2)(r - r0)^2. The Fresnel   *
         *  approximation is psi = (pi/2) * (r - r0)^2 / F^2.                 */
        const double psi = psi_factor * x_arr[m];

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ipsi = tmpl_CDouble_Polar(w_func[m], -psi);

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. That is,    *
         *  if rho = rho0 + x, then (rho - rho0)^2 = x^2. Similarly, if       *
         *  rho = rho0 - x, then (rho - rho0)^2 = x^2. Thus rho0 + x and      *
         *  rho0 - x will produce the same psi value. The Riemann sum then    *
         *  factors as:                                                       *
         *                                                                    *
         *        T_in[left] exp(-i psi) + T_in[right] exp(-i psi)            *
         *      = (T_in[left] + T_in[right]) exp(-i psi)                      *
         *                                                                    *
         *  which saves us some computational time. First, compute the        *
         *  T_in[left] + T_in[right] expression.                              */
        integrand = tmpl_CDouble_Add(tau->T_in[center-n], tau->T_in[center+n]);

        /*  Scaling by exp(-i psi) gives us the integrand for the transform.  *
         *  Add this to the partial sum (T_out) to compute the integral.      */
        tmpl_CDouble_MultiplyBy(&integrand, &exp_negative_ipsi);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        /*  Decrement the offset index to the next point. We started at the   *
         *  edge of the window and are moving inwards towards the center.     */
        n--;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    real = scale * (tau->T_out[center].dat[0] - tau->T_out[center].dat[1]);
    imag = scale * (tau->T_out[center].dat[0] + tau->T_out[center].dat[1]);
    tau->T_out[center].dat[0] = real;
    tau->T_out[center].dat[1] = imag;
}
/*  End of rssringoccs_Fresnel_Transform.                                     */
