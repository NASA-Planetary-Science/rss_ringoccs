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

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  size_t typedef given here.                                                */
#include <stddef.h>

extern void
rssringoccs_Fresnel_Transform(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Normalized(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
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

extern void
rssringoccs_Fresnel_Transform_Newton(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Norm(
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
rssringoccs_Fresnel_Transform_Perturbed_Newton_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Quartic(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Quartic_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Elliptical(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_Elliptical_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

#endif
/*  End of include guard.                                                     */
