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
 *      Typedef for the Tau object with geometry and diffraction data.        *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       December 28, 2025                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_TYPES_TAUOBJ_H
#define RSS_RINGOCCS_TYPES_TAUOBJ_H

/*  Booleans (True and False) are given here.                                 */
#include <libtmpl/include/tmpl_bool.h>

/*  Complex numbers found here.                                               */
#include <libtmpl/include/types/tmpl_complex_double.h>

/*  Typedef for window functions (function pointers f:R^2 -> R).              */
#include <libtmpl/include/types/tmpl_window_function_double.h>

/*  The psitype enum is defined here.                                         */
#include <rss_ringoccs/include/types/rss_ringoccs_psitype.h>

/*  DLP object typedef provided here.                                         */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  Structure that contains all of the necessary data. This includes geometry *
 *  data, diffraction data, diffraction corrected data, and forward modeling  *
 *  data.                                                                     */
typedef struct rssringoccs_TAUObj_Def {
    rssringoccs_DLPObj *dlp;
    tmpl_ComplexDouble *T_in;
    tmpl_ComplexDouble *T_out;
    tmpl_ComplexDouble *T_fwd;
    double *F_km_vals;
    double *k_vals;
    double *w_km_vals;
    double *tau_threshold_vals;
    double normeq;
    double sigma;
    double eccentricity;
    double periapse;
    double resolution_km;
    double perturb[5];
    double range[2];
    double root_finding_epsilon;
    unsigned int root_finding_max_iters;
    size_t start;
    size_t n_used;
    tmpl_WindowFunctionDouble window_func;
    enum rssringoccs_PsiType psinum;
    tmpl_Bool use_norm;
    tmpl_Bool use_fwd;
    tmpl_Bool bfac;
    tmpl_Bool verbose;
    tmpl_Bool error_occurred;
    const char *error_message;
    unsigned int order;
} rssringoccs_TAUObj;

#endif
/*  End of include guard.                                                     */
