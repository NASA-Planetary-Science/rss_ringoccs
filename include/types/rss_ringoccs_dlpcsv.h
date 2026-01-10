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
 *      Provides a struct for the data in a DLP.TAB file.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 5, 2026                                               *
 ******************************************************************************/

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_DLPCSV_H
#define RSS_RINGOCCS_TYPES_DLPCSV_H

/*  Booleans found here.                                                      */
#include <libtmpl/include/tmpl_bool.h>

/*  History object typedef is here. Each CSV object gets its own history.     */
#include <rss_ringoccs/include/types/rss_ringoccs_history.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Data structure for the DLP.TAB files on the PDS.                          */
typedef struct rssringoccs_DLPCSV_Def {
    double *rho_km_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *phi_rl_deg_vals;
    double *phi_ora_deg_vals;
    double *p_norm_vals;
    double *raw_tau_vals;
    double *phase_deg_vals;
    double *raw_tau_threshold_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *B_deg_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool use_deprecated;
    tmpl_Bool error_occurred;
    const char *error_message;
} rssringoccs_DLPCSV;

#endif
/*  End of include guard.                                                     */
