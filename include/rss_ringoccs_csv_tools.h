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

#include <rss_ringoccs/include/rss_ringoccs_history.h>

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
    double *obs_spacecraft_lat_deg_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool use_deprecated;
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
    rssringoccs_History *history;
    tmpl_Bool use_deprecated;
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
    rssringoccs_History *history;
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
    rssringoccs_History *history;
    tmpl_Bool use_deprecated;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_TauCSV;

/*  Data structure that contains all of the data from all four CSV formats    *
 *  interpolated so that the values are a function of radius, not time.       */
typedef struct rssringoccs_CSVData_Def {
    rssringoccs_GeoCSV *geo;
    rssringoccs_CalCSV *cal;
    rssringoccs_DLPCSV *dlp;
    rssringoccs_TauCSV *tau;
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
    double *tau_phase_deg_vals;
    double *tau_power_vals;
    double *tau_vals;
    size_t n_elements;
    rssringoccs_History *history;
    size_t geo_increment;
    size_t geo_decrement;
    tmpl_Bool use_deprecated;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_CSVData;

/*  Voyager / Uranus DLP files. Differ from the Cassini / Saturn ones.        */
typedef struct rssringoccs_UranusDLPCSV_Def {
    double *rho_km_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *phi_rl_deg_vals;
    double *phi_deg_vals;
    double *p_norm_vals;
    double *raw_tau_vals;
    double *phase_deg_vals;
    double *raw_tau_threshold_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *B_deg_vals;
    double *rho_dot_kms_vals;
    double *F_km_vals;
    double *D_km_vals;
    double *f_sky_hz_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool in_radians;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_UranusDLPCSV;

/*  Data structure that contains all of the data from all four CSV formats    *
 *  interpolated so that the values are a function of radius, not time.       */
typedef struct rssringoccs_UranusCSVData_Def {
    rssringoccs_GeoCSV *geo;
    rssringoccs_UranusDLPCSV *dlp;
    rssringoccs_TauCSV *tau;
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
    double *tau_phase_deg_vals;
    double *tau_power_vals;
    double *tau_vals;
    size_t n_elements;
    rssringoccs_History *history;
    size_t geo_increment;
    size_t geo_decrement;
    tmpl_Bool dlp_in_radians;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_UranusCSVData;

/*  Struct for the Merged CSV files (DLPM.TAB). This combines DLP, GEO, and   *
 *  CAL data, and can optionally contain TAU data as well.                    */
typedef struct rssringoccs_MergedCSVData_Def {
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
    double *tau_phase_deg_vals;
    double *tau_power_vals;
    double *tau_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_MergedCSVData;

/******************************************************************************
 *                               Cal CSV Tools                                *
 ******************************************************************************/
extern rssringoccs_CalCSV *rssringoccs_CalCSV_Extract(const char *filename);

extern void
rssringoccs_CalCSV_Check_Column_Count(rssringoccs_CalCSV *cal, FILE *fp);

extern void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV *cal);
extern void rssringoccs_CalCSV_Destroy(rssringoccs_CalCSV **cal);
extern void rssringoccs_CalCSV_Init(rssringoccs_CalCSV *cal);
extern void rssringoccs_CalCSV_Malloc(rssringoccs_CalCSV *cal, FILE *fp);
extern void rssringoccs_CalCSV_Read_Data(rssringoccs_CalCSV *cal, FILE *fp);

/******************************************************************************
 *                               DLP CSV Tools                                *
 ******************************************************************************/
extern rssringoccs_DLPCSV *
rssringoccs_DLPCSV_Extract(const char *filename, tmpl_Bool use_deprecated);

extern void
rssringoccs_DLPCSV_Check_Column_Count(rssringoccs_DLPCSV *dlp, FILE *fp);

extern void rssringoccs_DLPCSV_Destroy_Members(rssringoccs_DLPCSV *dlp);
extern void rssringoccs_DLPCSV_Destroy(rssringoccs_DLPCSV **dlp);
extern void rssringoccs_DLPCSV_Init(rssringoccs_DLPCSV *dlp);
extern void rssringoccs_DLPCSV_Malloc(rssringoccs_DLPCSV *dlp, FILE *fp);
extern void rssringoccs_DLPCSV_Read_Data(rssringoccs_DLPCSV *dlp, FILE *fp);

/******************************************************************************
 *                               Geo CSV Tools                                *
 ******************************************************************************/
extern rssringoccs_GeoCSV *
rssringoccs_GeoCSV_Extract(const char *filename, tmpl_Bool use_deprecated);

extern void
rssringoccs_GeoCSV_Check_Column_Count(rssringoccs_GeoCSV *geo, FILE *fp);

extern void rssringoccs_GeoCSV_Destroy_Members(rssringoccs_GeoCSV *geo);
extern void rssringoccs_GeoCSV_Destroy(rssringoccs_GeoCSV **geo);
extern void rssringoccs_GeoCSV_Init(rssringoccs_GeoCSV *geo);
extern void rssringoccs_GeoCSV_Malloc(rssringoccs_GeoCSV *geo, FILE *fp);
extern void rssringoccs_GeoCSV_Read_Data(rssringoccs_GeoCSV *geo, FILE *fp);

extern void
rssringoccs_GeoCSV_Write_History(rssringoccs_GeoCSV *geo, const char *filename);

/******************************************************************************
 *                               Tau CSV Tools                                *
 ******************************************************************************/
extern rssringoccs_TauCSV *
rssringoccs_TauCSV_Extract(const char *filename, tmpl_Bool use_deprecated);

extern void
rssringoccs_TauCSV_Check_Column_Count(rssringoccs_TauCSV *geo, FILE *fp);

extern void rssringoccs_TauCSV_Destroy_Members(rssringoccs_TauCSV *geo);
extern void rssringoccs_TauCSV_Destroy(rssringoccs_TauCSV **geo);
extern void rssringoccs_TauCSV_Init(rssringoccs_TauCSV *geo);
extern void rssringoccs_TauCSV_Malloc(rssringoccs_TauCSV *geo, FILE *fp);
extern void rssringoccs_TauCSV_Read_Data(rssringoccs_TauCSV *geo, FILE *fp);

/******************************************************************************
 *                              Merged CSV Tools                              *
 ******************************************************************************/
extern rssringoccs_MergedCSVData *
rssringoccs_MergedCSVData_Extract(const char *filename);

extern void
rssringoccs_MergedCSVData_Check_Column_Count(rssringoccs_MergedCSVData *dlpm,
                                             FILE *fp);

extern void
rssringoccs_MergedCSVData_Destroy_Members(rssringoccs_MergedCSVData *dlp);

extern void
rssringoccs_MergedCSVData_Destroy(rssringoccs_MergedCSVData **dlp);

extern void
rssringoccs_MergedCSVData_Init(rssringoccs_MergedCSVData *dlp);

extern void
rssringoccs_MergedCSVData_Malloc(rssringoccs_MergedCSVData *dlp, FILE *fp);

extern void
rssringoccs_MergedCSVData_Read_Data(rssringoccs_MergedCSVData *dlp, FILE *fp);

/******************************************************************************
 *                               CSV Data Tools                               *
 ******************************************************************************/
extern void
rssringoccs_CSVData_Extract_Geo(rssringoccs_CSVData *csv, const char *geo_file);

extern void
rssringoccs_CSVData_Extract_Cal(rssringoccs_CSVData *csv, const char *cal_file);

extern void
rssringoccs_CSVData_Extract_DLP(rssringoccs_CSVData *csv, const char *dlp_file);

extern void
rssringoccs_CSVData_Extract_Tau(rssringoccs_CSVData *csv, const char *tau_file);

extern void rssringoccs_CSVData_Init(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Destroy_Members(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Destroy(rssringoccs_CSVData **csv);
extern void rssringoccs_CSVData_Malloc(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Steal_DLP_Data(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Reverse_Geo_Variables(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Interpolate_Geo(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Interpolate_Cal(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Interpolate_Tau(rssringoccs_CSVData *csv);
extern void rssringoccs_CSVData_Check_Chord_Occ(rssringoccs_CSVData *csv);

extern rssringoccs_CSVData *
rssringoccs_CSVData_Extract(const char *geo,
                            const char *cal,
                            const char *dlp,
                            const char *tau,
                            tmpl_Bool use_deprecated);

/******************************************************************************
 *                            Uranus DLP CSV Tools                            *
 ******************************************************************************/
extern rssringoccs_UranusDLPCSV *
rssringoccs_UranusDLPCSV_Extract(const char *filename, tmpl_Bool in_radians);

extern void
rssringoccs_UranusDLPCSV_Check_Column_Count(rssringoccs_UranusDLPCSV *dlp,
                                            FILE *fp);

extern void
rssringoccs_UranusDLPCSV_Destroy_Members(rssringoccs_UranusDLPCSV *dlp);

extern void
rssringoccs_UranusDLPCSV_Destroy(rssringoccs_UranusDLPCSV **dlp);

extern void
rssringoccs_UranusDLPCSV_Init(rssringoccs_UranusDLPCSV *dlp);

extern void
rssringoccs_UranusDLPCSV_Malloc(rssringoccs_UranusDLPCSV *dlp, FILE *fp);

extern void
rssringoccs_UranusDLPCSV_Read_Data(rssringoccs_UranusDLPCSV *dlp, FILE *fp);

/******************************************************************************
 *                           Uranus CSV Data Tools                            *
 ******************************************************************************/
extern void
rssringoccs_UranusCSVData_Extract_Geo(rssringoccs_UranusCSVData *csv,
                                      const char *geo_file);

extern void
rssringoccs_UranusCSVData_Extract_DLP(rssringoccs_UranusCSVData *csv,
                                      const char *dlp_file);

extern void
rssringoccs_UranusCSVData_Extract_Tau(rssringoccs_UranusCSVData *csv,
                                      const char *tau_file);

extern void
rssringoccs_UranusCSVData_Init(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Destroy_Members(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Destroy(rssringoccs_UranusCSVData **csv);

extern void
rssringoccs_UranusCSVData_Malloc(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Steal_DLP_Data(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Reverse_Geo_Variables(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Interpolate_Geo(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Interpolate_Tau(rssringoccs_UranusCSVData *csv);

extern void
rssringoccs_UranusCSVData_Check_Chord_Occ(rssringoccs_UranusCSVData *csv);

extern rssringoccs_UranusCSVData *
rssringoccs_UranusCSVData_Extract(const char *geo,
                                  const char *dlp,
                                  const char *tau,
                                  tmpl_Bool dlp_in_radians);

#endif
/*  End of include guard.                                                     */
