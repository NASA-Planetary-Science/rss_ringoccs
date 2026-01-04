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

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_TAUOBJ_H
#define RSS_RINGOCCS_TYPES_TAUOBJ_H

/*  Typedef for window functions (function pointers f:R^2 -> R).              */
#include <libtmpl/include/types/tmpl_window_function_double.h>

/*  The psitype enum is defined here.                                         */
#include <rss_ringoccs/include/types/rss_ringoccs_psitype.h>

/*  Structure that contains all of the necessary data. This includes geometry *
 *  data, diffraction data, diffraction corrected data, and forward modeling  *
 *  data.                                                                     */
typedef struct rssringoccs_TAUObj_Def {
    tmpl_ComplexDouble *T_in;
    tmpl_ComplexDouble *T_out;
    tmpl_ComplexDouble *T_fwd;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_deg_vals;
    double *k_vals;
    double *rho_dot_kms_vals;
    double *B_deg_vals;
    double *D_km_vals;
    double *w_km_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *tau_threshold_vals;
    double *phi_rl_deg_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double dx_km;
    double normeq;
    double sigma;
    double eccentricity;
    double periapse;
    double resolution_km;
    double perturb[5];
    double rng_list[2];
    double rng_req[2];
    double EPS;
    unsigned int toler;
    size_t start;
    size_t n_used;
    size_t arr_size;
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

