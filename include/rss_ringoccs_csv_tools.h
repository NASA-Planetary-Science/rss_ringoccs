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
 *      Provide tools for extracting data from the .TAB files on the NASA PDS.*
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 30, 2020                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_CSV_TOOLS_H
#define RSS_RINGOCCS_CSV_TOOLS_H

/*  Boolean data types defined here.                                          */
#include <libtmpl/include/tmpl_bool.h>

/*  size_t typedef is given here.                                             */
#include <stddef.h>
#include <stdio.h>

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
    double *obs_spacecract_lat_deg_vals;
    size_t n_elements;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_GeoCSV;

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
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_DLPCSV;

/*  Data structure for the CAL.TAB files on the PDS.                          */
typedef struct rssringoccs_CalCSV_Def {
    double *t_oet_spm_vals;
    double *f_sky_pred_vals;
    double *f_sky_resid_fit_vals;
    double *p_free_vals;
    size_t n_elements;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_CalCSV;

/*  Data structure for the TAU.TAB files on the PDS.                          */
typedef struct rssringoccs_TauCSV_Def {
    double *rho_km_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *phi_rl_deg_vals;
    double *phi_ora_deg_vals;
    double *power_vals;
    double *tau_vals;
    double *phase_deg_vals;
    double *tau_threshold_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *B_deg_vals;
    size_t n_elements;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_TauCSV;

/*  Data structure that contains all of the data from all four CSV formats    *
 *  interpolated so that the values are a function of radius, not time.       */
typedef struct rssringoccs_CSVData_Def {
    double *B_deg_vals;
    double *D_km_vals;
    double *f_sky_hz_vals;
    double *p_norm_vals;
    double *raw_tau_vals;
    double *phase_deg_vals;
    double *phi_deg_vals;
    double *phi_rl_deg_vals;
    double *raw_tau_threshold_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *rho_dot_kms_vals;
    double *rho_km_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *tau_phase;
    double *tau_power;
    double *tau_vals;
    size_t n_elements;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_CSVData;

extern rssringoccs_GeoCSV *
rssringoccs_Get_Geo(const char *filename, tmpl_Bool use_deprecated);

extern void rssringoccs_Destroy_GeoCSV_Members(rssringoccs_GeoCSV *geo);
extern void rssringoccs_Destroy_GeoCSV(rssringoccs_GeoCSV **geo);

extern rssringoccs_DLPCSV *
rssringoccs_Get_DLP(const char *filename, tmpl_Bool use_deprecated);

extern void rssringoccs_Destroy_DLPCSV_Members(rssringoccs_DLPCSV *dlp);
extern void rssringoccs_Destroy_DLPCSV(rssringoccs_DLPCSV **dlp);

extern rssringoccs_CalCSV *rssringoccs_Get_Cal(const char *filename);
extern void rssringoccs_Destroy_CalCSV_Members(rssringoccs_CalCSV *cal);
extern void rssringoccs_Destroy_CalCSV(rssringoccs_CalCSV **cal);

extern rssringoccs_TauCSV *
rssringoccs_Get_Tau(const char *filename, tmpl_Bool use_deprecated);

extern void rssringoccs_Destroy_TauCSV_Members(rssringoccs_TauCSV *dlp);
extern void rssringoccs_Destroy_TauCSV(rssringoccs_TauCSV **dlp);

extern rssringoccs_CSVData *
rssringoccs_Extract_CSV_Data(const char *geo,
                             const char *cal,
                             const char *dlp,
                             const char *tau,
                             tmpl_Bool use_deprecated);

extern void rssringoccs_Destroy_CSV_Members(rssringoccs_CSVData *csv);

extern void rssringoccs_Destroy_CSV(rssringoccs_CSVData **csv);

#endif
/*  End of include guard.                                                     */
