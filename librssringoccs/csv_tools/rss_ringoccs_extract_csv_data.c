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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <rss_ringoccs/include/rss_ringoccs_interpolate.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define __MALLOC_CSV_VAR__(var)                                                \
    csv_data->var = malloc(sizeof(*csv_data->var)*csv_data->n_elements);       \
    if (csv_data->var == NULL)                                                 \
    {                                                                          \
        csv_data->error_occurred = rssringoccs_True;                           \
        csv_data->error_message = rssringoccs_strdup(                          \
            "Error Encountered: rss_ringoccs\n"                                \
            "\trssringoccs_Extract_CSV_Data\n\n"                               \
            "Malloc returned NULL for csv_data member. Aborting.\n"            \
        );                                                                     \
        rssringoccs_Destroy_GeoCSV(&geo_dat);                                  \
        rssringoccs_Destroy_DLPCSV(&dlp_dat);                                  \
        rssringoccs_Destroy_CalCSV(&cal_dat);                                  \
        rssringoccs_Destroy_TauCSV(&tau_dat);                                  \
        rssringoccs_Destroy_CSV_Members(csv_data);                             \
        return csv_data;                                                       \
    }

rssringoccs_CSVData *
rssringoccs_Extract_CSV_Data(const char *geo,
                             const char *cal,
                             const char *dlp,
                             const char *tau,
                             rssringoccs_Bool use_deprecated)
{
    rssringoccs_GeoCSV  *geo_dat;
    rssringoccs_DLPCSV  *dlp_dat;
    rssringoccs_CalCSV  *cal_dat;
    rssringoccs_TauCSV  *tau_dat;
    rssringoccs_CSVData *csv_data;
    unsigned long n;
    double min_dr_dt, max_dr_dt, temp;
    double *geo_rho, *geo_rho_dot, *geo_D;

    csv_data = malloc(sizeof(*csv_data));

    if (csv_data == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Extract_CSV_Data\n\n"
             "Malloc failed and returned NULL for csv_data. Returning.\n");
        return NULL;
    }

    /*  Initialize the members to NULL. This will prevent functions from      *
     *  trying to free pointers that weren't malloc'd in the event of error.  */
    csv_data->B_rad_vals = NULL;
    csv_data->D_km_vals = NULL;
    csv_data->f_sky_hz_vals = NULL;
    csv_data->p_norm_vals = NULL;
    csv_data->power_vals = NULL;
    csv_data->phase_rad_vals = NULL;
    csv_data->phase_vals = NULL;
    csv_data->phi_rad_vals = NULL;
    csv_data->phi_rl_rad_vals = NULL;
    csv_data->raw_tau_threshold_vals = NULL;
    csv_data->rho_corr_pole_km_vals = NULL;
    csv_data->rho_corr_timing_km_vals = NULL;
    csv_data->rho_dot_kms_vals = NULL;
    csv_data->rho_km_vals = NULL;
    csv_data->rx_km_vals = NULL;
    csv_data->ry_km_vals = NULL;
    csv_data->rz_km_vals = NULL;
    csv_data->t_oet_spm_vals = NULL;
    csv_data->t_ret_spm_vals = NULL;
    csv_data->t_set_spm_vals = NULL;
    csv_data->tau_rho = NULL;
    csv_data->tau_power = NULL;
    csv_data->tau_vals = NULL;
    csv_data->error_message = NULL;

    geo_dat = rssringoccs_Get_Geo(geo, use_deprecated);
    if (geo_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned NULL for geo_dat. Aborting.\n"
        );
        return csv_data;
    }
    else if (geo_dat->n_elements == 0)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }

    dlp_dat = rssringoccs_Get_DLP(dlp, use_deprecated);
    if (dlp_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned NULL for dlp_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }
    else if (dlp_dat->n_elements == 0)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }

    cal_dat = rssringoccs_Get_Cal(cal);
    if (cal_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned NULL for cal_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }
    else if (cal_dat->n_elements == 0)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }

    tau_dat = rssringoccs_Get_Tau(tau, use_deprecated);
    if (tau_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Tau returned NULL for tau_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }
    else if (tau_dat->n_elements == 0)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_TAU returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        rssringoccs_Destroy_TauCSV(&tau_dat);
        return csv_data;
    }

    /*  Grab the number of elements from the DLP CSV. This will be the number *
     *  of elements in the output.                                            */
    csv_data->n_elements = dlp_dat->n_elements;

    /*  The __MALLOC_CSV_VAR__ macro contains an if-then statement with       *
     *  braces {} hence there is no need for a semi-colon at the end.         */
    __MALLOC_CSV_VAR__(rho_km_vals)
    __MALLOC_CSV_VAR__(raw_tau_vals)
    __MALLOC_CSV_VAR__(phase_rad_vals)
    __MALLOC_CSV_VAR__(phi_rad_vals)
    __MALLOC_CSV_VAR__(B_rad_vals)
    __MALLOC_CSV_VAR__(t_oet_spm_vals)
    __MALLOC_CSV_VAR__(t_ret_spm_vals)
    __MALLOC_CSV_VAR__(t_set_spm_vals)
    __MALLOC_CSV_VAR__(rho_corr_pole_km_vals)
    __MALLOC_CSV_VAR__(rho_corr_timing_km_vals)
    __MALLOC_CSV_VAR__(phi_rl_rad_vals)
    __MALLOC_CSV_VAR__(raw_tau_threshold_vals)

    for (n=0; n<csv_data->n_elements; ++n)
    {
        csv_data->rho_km_vals[n]  = dlp_dat->rho_km_vals[n];
        csv_data->raw_tau_vals[n] = dlp_dat->raw_tau_vals[n];
        csv_data->phase_rad_vals[n]
            = rssringoccs_Deg_To_Rad*dlp_dat->phase_deg_vals[n];
        csv_data->phi_rad_vals[n]
            = rssringoccs_Deg_To_Rad*dlp_dat->phi_ora_deg_vals[n];
        csv_data->B_rad_vals[n]
            = rssringoccs_Deg_To_Rad*dlp_dat->B_deg_vals[n];
        csv_data->t_ret_spm_vals[n] = dlp_dat->t_ret_spm_vals[n];
        csv_data->t_set_spm_vals[n] = dlp_dat->t_set_spm_vals[n];
        csv_data->t_oet_spm_vals[n] = dlp_dat->t_oet_spm_vals[n];
        csv_data->rho_corr_pole_km_vals[n]
            = dlp_dat->rho_corr_pole_km_vals[n];
        csv_data->rho_corr_timing_km_vals[n]
            = dlp_dat->rho_corr_timing_km_vals[n];
        csv_data->phi_rl_rad_vals[n]
            = rssringoccs_Deg_To_Rad*dlp_dat->phi_rl_deg_vals[n];
        csv_data->raw_tau_threshold_vals[n]
            = dlp_dat->raw_tau_threshold_vals[n];
    }

    geo_rho     = geo_dat->rho_km_vals;
    geo_D       = geo_dat->D_km_vals;
    geo_rho_dot = geo_dat->rho_dot_kms_vals;

    min_dr_dt = -rssringoccs_Infinity;
    max_dr_dt = rssringoccs_Infinity;
    for (n=0; n<csv_data->n_elements-1; ++n)
    {
        temp = (csv_data->rho_km_vals[n+1] - csv_data->rho_km_vals[n])   /
               (csv_data->t_set_spm_vals[n+1] - csv_data->t_set_spm_vals[n]);
        if (temp < min_dr_dt)
            min_dr_dt = temp;

        if (max_dr_dt < temp)
            max_dr_dt = temp;
    }

    if ((min_dr_dt < 0.0) && (max_dr_dt > 0.0))
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "\rdrho/dt has positive and negative values. Check your DLP file.\n"
            "\rIt is likely a chord occultation and needs to be split into\n"
            "\ringress and egress portions.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        rssringoccs_Destroy_TauCSV(&tau_dat);
        rssringoccs_Destroy_CSV_Members(csv_data);
    }
    else if ((min_dr_dt == 0.0) || (max_dr_dt == 0.0))
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "\rdrho/dt has zero-valued elements. Check your DLP file.\n"
            "\rIt is likely a chord occultation and needs to be split into\n"
            "\ringress and egress portions.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        rssringoccs_Destroy_TauCSV(&tau_dat);
        rssringoccs_Destroy_CSV_Members(csv_data);
    }
    else if (max_dr_dt < 0.0)
    {
        rssringoccs_Reverse_Double_Array(geo_rho, geo_dat->n_elements);
        rssringoccs_Reverse_Double_Array(geo_rho_dot, geo_dat->n_elements);
        rssringoccs_Reverse_Double_Array(geo_D, geo_dat->n_elements);
    }

    rssringoccs_Destroy_GeoCSV(&geo_dat);
    rssringoccs_Destroy_DLPCSV(&dlp_dat);
    rssringoccs_Destroy_CalCSV(&cal_dat);
    rssringoccs_Destroy_TauCSV(&tau_dat);

    return csv_data;

}
