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
 *      Provides a struct for the data in a GEO.TAB file.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 5, 2026                                               *
 ******************************************************************************/

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_GEOCSV_H
#define RSS_RINGOCCS_TYPES_GEOCSV_H

/*  Booleans found here.                                                      */
#include <libtmpl/include/tmpl_bool.h>

/*  History object typedef is here. Each CSV object gets its own history.     */
#include <rss_ringoccs/include/types/rss_ringoccs_history.h>

/*  Data structure for the GEO.TAB files on the PDS.                          */
typedef struct rssringoccs_GeoCSV_Def {
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_km_vals;
    double *phi_rl_deg_vals;
    double *phi_ora_deg_vals;
    double *B_deg_vals;
    double *D_km_vals;
    double *rho_dot_kms_vals;
    double *phi_rl_dot_kms_vals;
    double *F_km_vals;
    double *R_imp_km_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double *vx_kms_vals;
    double *vy_kms_vals;
    double *vz_kms_vals;
    double *obs_spacecraft_lat_deg_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool use_deprecated;
    tmpl_Bool error_occurred;
    const char *error_message;
} rssringoccs_GeoCSV;

#endif
/*  End of include guard.                                                     */
