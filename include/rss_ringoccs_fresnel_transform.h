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
rssringoccs_Fresnel_Transform_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Legendre_Even(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    const double * TMPL_RESTRICT const coeffs,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Legendre_Even_Norm(
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
rssringoccs_Fresnel_Transform_Newton_D(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Old(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Old_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_dD_dphi(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_dD_dphi_Norm(
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

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Transform_Newton_Quadratic                        *
 *  Purpose:                                                                  *
 *      Performs the Fresnel inverse integral by approximating the stationary *
 *      Fresnel kernel with a quadratic polynomial.                           *
 *  Arguments:                                                                *
 *      rssringoccs_TAUObj *tau:                                              *
 *          A pointer to the tau object with geometry, power, and phase data. *
 *          The reconstructed profie is stored here as well.                  *
 *      const double *w_func:                                                 *
 *          An array corresponding to the values of the window function       *
 *          across the data set.                                              *
 *      size_t n_pts:                                                         *
 *          The number of points in the window.                               *
 *      size_t center:                                                        *
 *          The index corresponding to the center of the window.              *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Source:                                                                   *
 *      src/fresnel_transform/rss_ringoccs_fresnel_transform_quadratic.c      *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Transform_Newton_Quadratic(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Fresnel_Transform_Newton_Quadratic_Norm                   *
 *  Purpose:                                                                  *
 *      Same as rssringoccs_Fresnel_Transform_Quadratic but normalizes the    *
 *      output by the size of the window.                                     *
 *  Arguments:                                                                *
 *      rssringoccs_TAUObj *tau:                                              *
 *          A pointer to the tau object with geometry, power, and phase data. *
 *          The reconstructed profie is stored here as well.                  *
 *      const double *w_func:                                                 *
 *          An array corresponding to the values of the window function       *
 *          across the data set.                                              *
 *      size_t n_pts:                                                         *
 *          The number of points in the window.                               *
 *      size_t center:                                                        *
 *          The index corresponding to the center of the window.              *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Source:                                                                   *
 *      src/fresnel_transform/rss_ringoccs_fresnel_transform_quadratic_norm.c *
 ******************************************************************************/
extern void
rssringoccs_Fresnel_Transform_Newton_Quadratic_Norm(
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
rssringoccs_Fresnel_Transform_Newton_D_Quartic(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Quartic_Norm(
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
