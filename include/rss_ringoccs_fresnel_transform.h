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
 *  Purpose:                                                                  *
 *      Provides various methods of Fresnel inversion for radio occultations. *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       December 15, 2020                                             *
 ******************************************************************************/

/*  Include guard to avoid including this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_TRANFORM_H
#define RSS_RINGOCCS_FRESNEL_TRANFORM_H

#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/types/tmpl_complex_double.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef given here.                                                */
#include <stddef.h>

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Kernel                                            *
 *  Purpose:                                                                  *
 *      Computes the Fresnel kernel from the geometry data in a Tau object.   *
 *  Arguments:                                                                *
 *      tau (const rssringoccs_TAUObj * const):                               *
 *          The Tau object with the geometry and diffraction data.            *
 *      center (const size_t):                                                *
 *          The index for the center of the window of integration.            *
 *      offset (const size_t):                                                *
 *          The index for the dummy variable of integration.                  *
 *  Output:                                                                   *
 *      kernel (tmpl_ComplexDouble).                                          *
 *  Notes:                                                                    *
 *      1.) There are no checks for NULL pointers. It is assumed tau is valid.*
 *                                                                            *
 *      2.) There are no checks for errors. This function is called inside    *
 *          the main for-loop of the various transforms. It is the users      *
 *          obligation to ensure the input data is valid.                     *
 *                                                                            *
 *      3.) There are no checks that the center and offset arguments          *
 *          correspond to indices for the actual data.                        *
 *                                                                            *
 *      4.) It is assumed that the input Tau object has all of its angles in  *
 *          degrees. This should always be the case, regardless.              *
 ******************************************************************************/
extern tmpl_ComplexDouble
rssringoccs_Fresnel_Kernel(const rssringoccs_TAUObj * const tau,
                           const size_t center,
                           const size_t offset);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Phase_And_Weight                                  *
 *  Purpose:                                                                  *
 *      Computes the Fresnel kernel from the geometry data in a Tau object.   *
 *      The output is in polar form (magnitude, angle).                       *
 *  Arguments:                                                                *
 *      tau (const rssringoccs_TAUObj * const):                               *
 *          The Tau object with the geometry and diffraction data.            *
 *      center (const size_t):                                                *
 *          The index for the center of the window of integration.            *
 *      offset (const size_t):                                                *
 *          The index for the dummy variable of integration.                  *
 *      weight (double * TMPL_RESTRICT const weight):                         *
 *          The scale factor for the Fresnel kernel.                          *
 *      psi (double * TMPL_RESTRICT const):                                   *
 *          The Fresnel phase, the part that is complex exponentiated.        *
 *          This is in radians.                                               *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      1.) There are no checks for NULL pointers. It is assumed tau is valid.*
 *                                                                            *
 *      2.) There are no checks for errors. This function is called inside    *
 *          the main for-loop of the various transforms. It is the users      *
 *          obligation to ensure the input data is valid.                     *
 *                                                                            *
 *      3.) There are no checks that the center and offset arguments          *
 *          correspond to indices for the actual data.                        *
 *                                                                            *
 *      4.) It is assumed that the input Tau object has all of its angles in  *
 *          degrees. This should always be the case, regardless.              *
 *                                                                            *
 *      5.) Note, the output angle is in radians, not degrees. This is        *
 *          because the Newton-Raphson method used to compute the Fresnel     *
 *          phase works in radians (implemented in libtmpl).                  *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Phase_And_Weight(
    const rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const size_t center,
    const size_t offset,
    double * TMPL_RESTRICT const weight,
    double * TMPL_RESTRICT const psi
);

extern void
rssringoccs_Fresnel_Transform(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const size_t center,
    const size_t n_pts
);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Transform_Even_Polynomial                         *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using an even degree polynomial       *
 *      approximation for the stationary Fresnel kernel.                      *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          A pointer to a Tau object. This contains all of the geometery and *
 *          diffraction data. This function writes the newly corrected data   *
 *          to tau->T_out[center] (see below for the description of center).  *
 *      x_arr (const double * TMPL_RESTRICT const):                           *
 *          The array r[n] - r[center] as n varies between 0 and n_pts, where *
 *          r is the radius. The array thus starts at the left-most endpoint  *
 *          of the window and varies to the center. x_arr[n_pts + 1] is the   *
 *          center of the window.                                             *
 *      w_funct (const double * TMPL_RESTRICT const):                         *
 *          The pre-computed window function as a function of x_arr. That is, *
 *          w_func[n] is the value of the window function at x_arr[n].        *
 *      coeffs (const double * TMPL_RESTRICT const):                          *
 *          The coefficients for the polynomial approximation. There must be  *
 *          at least tau->order elements to the array. coeffs[0] represents   *
 *          the constant coefficients, coeffs[tau->order - 1] corresponds to  *
 *          the coefficients of the highest order term.                       *
 *      n_pts (size_t):                                                       *
 *          The number of points in the x_arr array. There are 2 * n_pts + 1  *
 *          points in the window, n_pts to the left of center, n_pts to the   *
 *          right, and one point in the center.                               *
 *      center (size_t):                                                      *
 *          The index corresponding to the center of the window. There must   *
 *          be n_pts to the left and right of center in the data.             *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      1.) Since the coefficient array is pre-computed, you may use any      *
 *          polynomial approximation you wish. rss_ringoccs mostly uses this  *
 *          function with Legendre polynomials.                               *
 *      2.) Various polynomial approximations make certain assumptions. The   *
 *          Legendre polynomials assume the first iteration of Newton's       *
 *          method produces a decent approximation for the stationary         *
 *          azimuth angle. This is not true for low ring openging angles.     *
 *      3.) There are no checks for NULL pointers. You are responsible for    *
 *          using this function with valid pointers.                          *
 *      4.) This function does not normalize by the window width. Avoid using *
 *          with coarse resolutions.                                          *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Transform_Even_Polynomial(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Transform_Normalized_Even_Polynomial              *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using an even degree polynomial       *
 *      approximation for the stationary Fresnel kernel, normalizing the      *
 *      end result by the window width.                                       *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          A pointer to a Tau object. This contains all of the geometery and *
 *          diffraction data. This function writes the newly corrected data   *
 *          to tau->T_out[center] (see below for the description of center).  *
 *      x_arr (const double * TMPL_RESTRICT const):                           *
 *          The array r[n] - r[center] as n varies between 0 and n_pts, where *
 *          r is the radius. The array thus starts at the left-most endpoint  *
 *          of the window and varies to the center. x_arr[n_pts + 1] is the   *
 *          center of the window.                                             *
 *      w_funct (const double * TMPL_RESTRICT const):                         *
 *          The pre-computed window function as a function of x_arr. That is, *
 *          w_func[n] is the value of the window function at x_arr[n].        *
 *      coeffs (const double * TMPL_RESTRICT const):                          *
 *          The coefficients for the polynomial approximation. There must be  *
 *          at least tau->order elements to the array. coeffs[0] represents   *
 *          the constant coefficients, coeffs[tau->order - 1] corresponds to  *
 *          the coefficients of the highest order term.                       *
 *      n_pts (size_t):                                                       *
 *          The number of points in the x_arr array. There are 2 * n_pts + 1  *
 *          points in the window, n_pts to the left of center, n_pts to the   *
 *          right, and one point in the center.                               *
 *      center (size_t):                                                      *
 *          The index corresponding to the center of the window. There must   *
 *          be n_pts to the left and right of center in the data.             *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      1.) Since the coefficient array is pre-computed, you may use any      *
 *          polynomial approximation you wish. rss_ringoccs mostly uses this  *
 *          function with Legendre polynomials.                               *
 *      2.) Various polynomial approximations make certain assumptions. The   *
 *          Legendre polynomials assume the first iteration of Newton's       *
 *          method produces a decent approximation for the stationary         *
 *          azimuth angle. This is not true for low ring openging angles.     *
 *      3.) There are no checks for NULL pointers. You are responsible for    *
 *          using this function with valid pointers.                          *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Transform_Normalized_Even_Polynomial(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Legendre_Odd(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Legendre_Odd_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Transform_Newton_Riemann                          *
 *  Purpose:                                                                  *
 *      Performs diffraction reconstruction using Newton's method with a left *
 *      Riemann sum for numerical integration.                                *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          A pointer to a Tau object. This contains all of the geometery and *
 *          diffraction data. This function writes the newly corrected data   *
 *          to tau->T_out[center] (see below for the description of center).  *
 *      w_funct (const double * TMPL_RESTRICT const):                         *
 *          The pre-computed window function as a function of x_arr. That is, *
 *          w_func[n] is the value of the window function at x_arr[n].        *
 *      n_pts (size_t):                                                       *
 *          The number of points in the window.                               *
 *      center (size_t):                                                      *
 *          The index corresponding to the center of the window. There must   *
 *          be n_pts to the left and right of center in the data.             *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      1.) The Riemann sum produces poor results when the Fresnel phase has  *
 *          a large derivative with respect to rho. That is, when dpsi / drho *
 *          is large (say, larger than 2 cycles per bin), this method         *
 *          produces a poor reconstruction. Use the Filon methods instead.    *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Transform_Newton_Riemann(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Filon01(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Filon11(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Filon02(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Elliptical_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Perturbed_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton4(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton8(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton16(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

#endif
/*  End of include guard.                                                     */
