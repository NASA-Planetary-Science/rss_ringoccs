/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                         rss_ringoccs_csv_tools                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide tools for extracting data from the .TAB files on the NASA PDS.*
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 30, 2020                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_CSV_TOOLS_H__
#define __RSS_RINGOCCS_CSV_TOOLS_H__

/*  Boolean data types defined here.                                          */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

/*  Data structure for the GEO.TAB files on the PDS.                          */
typedef struct rssringoccs_GeoCSV {
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
    double *obs_spacecract_lat_deg_vals;
    unsigned long n_elements;
    rssringoccs_Bool error_occurred;
    char *error_message;
} rssringoccs_GeoCSV;

extern rssringoccs_GeoCSV *
rssringoccs_Get_Geo(const char *filename, rssringoccs_Bool use_deprecated);

extern void rssringoccs_Destroy_GeoCSV(rssringoccs_GeoCSV **geo);

#endif
/*  End of include guard.                                                     */
