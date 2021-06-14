/******************************************************************************
 *                                 LICENSE                                    *
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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_interpolate.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>
#include <math.h>

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_CSV_VAR(var)                                                    \
    csv_data->var = malloc(sizeof(*csv_data->var)*csv_data->n_elements);       \
    if (csv_data->var == NULL)                                                 \
    {                                                                          \
        csv_data->error_occurred = tmpl_True;                                  \
        csv_data->error_message = tmpl_strdup(                                 \
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
                             tmpl_Bool use_deprecated)
{
    /*  Pointers for the Geo, Cal, DLP, and Tau CSV objects.                  */
    rssringoccs_GeoCSV  *geo_dat;
    rssringoccs_DLPCSV  *dlp_dat;
    rssringoccs_CalCSV  *cal_dat;
    rssringoccs_TauCSV  *tau_dat;

    /*  A pointer to the CSV object.                                          */
    rssringoccs_CSVData *csv_data;

    /*  Variable for indexing.                                                */
    unsigned long n;

    /*  Variables for checking the geometry of the occultation.               */
    double min_dr_dt, max_dr_dt, temp, raw_mu;
    double *cal_f_sky_hz_vals;

    /*  Allocate memory for the CSV object.                                   */
    csv_data = malloc(sizeof(*csv_data));

    /*  Check if malloc failed.                                               */
    if (csv_data == NULL)
        return NULL;

    /*  Initialize the members to NULL. This will prevent functions from      *
     *  trying to free pointers that weren't malloc'd in the event of error.  */
    csv_data->B_rad_vals = NULL;
    csv_data->D_km_vals = NULL;
    csv_data->f_sky_hz_vals = NULL;
    csv_data->p_norm_vals = NULL;
    csv_data->raw_tau_vals = NULL;
    csv_data->phase_rad_vals = NULL;
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
    csv_data->tau_phase = NULL;
    csv_data->tau_vals = NULL;
    csv_data->error_message = NULL;
    csv_data->error_occurred = tmpl_False;

    /*  Extract the data from the GEO.TAB file.                               */
    geo_dat = rssringoccs_Get_Geo(geo, use_deprecated);

    /*  Check for errors.                                                     */
    if (geo_dat == NULL)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned NULL for geo_dat. Aborting.\n"
        );
        return csv_data;
    }
    else if (geo_dat->n_elements == 0U)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }

    /*  Extract the data from the DLP.TAB file.                               */
    dlp_dat = rssringoccs_Get_DLP(dlp, use_deprecated);

    /*  Check for errors.                                                     */
    if (dlp_dat == NULL)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned NULL for dlp_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }
    else if (dlp_dat->n_elements == 0U)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }

    /*  Extract the data from the CAL.TAB file.                               */
    cal_dat = rssringoccs_Get_Cal(cal);

    /*  Check for errors.                                                     */
    if (cal_dat == NULL)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned NULL for cal_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }
    else if (cal_dat->n_elements == 0U)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned an empty struct. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }

    cal_f_sky_hz_vals = malloc(sizeof(*cal_f_sky_hz_vals)*cal_dat->n_elements);

    if (cal_f_sky_hz_vals == NULL)
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "malloc returned NULL for cal_f_sky_hz_vals. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }

    for (n = 0UL; n < cal_dat->n_elements; ++n)
        cal_f_sky_hz_vals[n] = cal_dat->f_sky_pred_vals[n] - 
                               cal_dat->f_sky_resid_fit_vals[n];

    /*  Extract the data from the TAU.TAB file.                               */
    if (tau)
    {
        tau_dat = rssringoccs_Get_Tau(tau, use_deprecated);

        /*  Check for errors.                                                     */
        if (tau_dat == NULL)
        {
            csv_data->error_occurred = tmpl_True;
            csv_data->error_message = tmpl_strdup(
                "Error Encountered: rss_ringoccs\n"
                "\trssringoccs_Extract_CSV_Data\n\n"
                "rssringoccs_Get_Tau returned NULL for tau_dat. Aborting.\n"
            );
            rssringoccs_Destroy_GeoCSV(&geo_dat);
            rssringoccs_Destroy_DLPCSV(&dlp_dat);
            rssringoccs_Destroy_CalCSV(&cal_dat);
            return csv_data;
        }
        else if (tau_dat->n_elements == 0U)
        {
            csv_data->error_occurred = tmpl_True;
            csv_data->error_message = tmpl_strdup(
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
    }
    else
        tau_dat = NULL;

    /*  Grab the number of elements from the DLP CSV. This will be the number *
     *  of elements in the output.                                            */
    csv_data->n_elements = dlp_dat->n_elements;

    /*  The MALLOC_CSV_VAR macro contains an if-then statement with           *
     *  braces {} hence there is no need for a semi-colon at the end.         */
    MALLOC_CSV_VAR(rho_km_vals)
    MALLOC_CSV_VAR(raw_tau_vals)
    MALLOC_CSV_VAR(p_norm_vals)
    MALLOC_CSV_VAR(phase_rad_vals)
    MALLOC_CSV_VAR(phi_rad_vals)
    MALLOC_CSV_VAR(B_rad_vals)
    MALLOC_CSV_VAR(t_oet_spm_vals)
    MALLOC_CSV_VAR(t_ret_spm_vals)
    MALLOC_CSV_VAR(t_set_spm_vals)
    MALLOC_CSV_VAR(rho_corr_pole_km_vals)
    MALLOC_CSV_VAR(rho_corr_timing_km_vals)
    MALLOC_CSV_VAR(phi_rl_rad_vals)
    MALLOC_CSV_VAR(raw_tau_threshold_vals)
    MALLOC_CSV_VAR(f_sky_hz_vals)
    MALLOC_CSV_VAR(D_km_vals)
    MALLOC_CSV_VAR(rho_dot_kms_vals)
    MALLOC_CSV_VAR(rx_km_vals)
    MALLOC_CSV_VAR(ry_km_vals)
    MALLOC_CSV_VAR(rz_km_vals)

    /*  Extract the DLP variables and store them in the CSV struct.           */
    for (n = 0U; n < csv_data->n_elements; ++n)
    {
        csv_data->rho_km_vals[n] = dlp_dat->rho_km_vals[n];
        csv_data->raw_tau_vals[n] = dlp_dat->raw_tau_vals[n];
        csv_data->phase_rad_vals[n]
            = tmpl_Deg_to_Rad*dlp_dat->phase_deg_vals[n];
        csv_data->phi_rad_vals[n]
            = tmpl_Deg_to_Rad*dlp_dat->phi_ora_deg_vals[n];
        csv_data->B_rad_vals[n] = tmpl_Deg_to_Rad*dlp_dat->B_deg_vals[n];
        csv_data->t_ret_spm_vals[n] = dlp_dat->t_ret_spm_vals[n];
        csv_data->t_set_spm_vals[n] = dlp_dat->t_set_spm_vals[n];
        csv_data->t_oet_spm_vals[n] = dlp_dat->t_oet_spm_vals[n];
        csv_data->rho_corr_pole_km_vals[n] = dlp_dat->rho_corr_pole_km_vals[n];
        csv_data->rho_corr_timing_km_vals[n]
            = dlp_dat->rho_corr_timing_km_vals[n];
        csv_data->phi_rl_rad_vals[n]
            = tmpl_Deg_to_Rad*dlp_dat->phi_rl_deg_vals[n];
        csv_data->raw_tau_threshold_vals[n]
            = dlp_dat->raw_tau_threshold_vals[n];
        raw_mu = sin(fabs(csv_data->B_rad_vals[n]));
        csv_data->p_norm_vals[n] = exp(-csv_data->raw_tau_vals[n] / raw_mu);
    }

    temp = (csv_data->rho_km_vals[1] - csv_data->rho_km_vals[0]) /
           (csv_data->t_set_spm_vals[1] - csv_data->t_set_spm_vals[0]);

    min_dr_dt = temp;
    max_dr_dt = temp;

    /*  Find the minimum and maximum of drho/dt.                              */
    for (n = 1U; n < csv_data->n_elements - 1U; ++n)
    {
        temp = (csv_data->rho_km_vals[n+1U] - csv_data->rho_km_vals[n])   /
               (csv_data->t_set_spm_vals[n+1U] - csv_data->t_set_spm_vals[n]);

        if (temp < min_dr_dt)
            min_dr_dt = temp;

        if (max_dr_dt < temp)
            max_dr_dt = temp;
    }

    /*  Check for errors.                                                     */
    if ((min_dr_dt < 0.0) && (max_dr_dt > 0.0))
    {
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
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
        csv_data->error_occurred = tmpl_True;
        csv_data->error_message = tmpl_strdup(
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
        tmpl_Reverse_Double_Array(geo_dat->rho_km_vals, geo_dat->n_elements);
        tmpl_Reverse_Double_Array(geo_dat->rho_dot_kms_vals,
                                  geo_dat->n_elements);
        tmpl_Reverse_Double_Array(geo_dat->D_km_vals, geo_dat->n_elements);
    }

    tmpl_Double_Sorted_Interp1d(geo_dat->rho_km_vals,
                                geo_dat->D_km_vals,
                                geo_dat->n_elements,
                                csv_data->rho_km_vals,
                                csv_data->D_km_vals,
                                csv_data->n_elements);

    tmpl_Double_Sorted_Interp1d(geo_dat->rho_km_vals,
                                geo_dat->rho_dot_kms_vals,
                                geo_dat->n_elements,
                                csv_data->rho_km_vals,
                                csv_data->rho_dot_kms_vals,
                                csv_data->n_elements);

    tmpl_Double_Sorted_Interp1d(geo_dat->rho_km_vals,
                                geo_dat->rx_km_vals,
                                geo_dat->n_elements,
                                csv_data->rho_km_vals,
                                csv_data->rx_km_vals,
                                csv_data->n_elements);

    tmpl_Double_Sorted_Interp1d(geo_dat->rho_km_vals,
                                geo_dat->ry_km_vals,
                                geo_dat->n_elements,
                                csv_data->rho_km_vals,
                                csv_data->ry_km_vals,
                                csv_data->n_elements);

    tmpl_Double_Sorted_Interp1d(geo_dat->rho_km_vals,
                                geo_dat->rz_km_vals,
                                geo_dat->n_elements,
                                csv_data->rho_km_vals,
                                csv_data->rz_km_vals,
                                csv_data->n_elements);

    /*  Free the Geo, Cal, DLP, and TAU structs.                              */
    rssringoccs_Destroy_GeoCSV(&geo_dat);
    rssringoccs_Destroy_DLPCSV(&dlp_dat);
    rssringoccs_Destroy_CalCSV(&cal_dat);
    rssringoccs_Destroy_TauCSV(&tau_dat);
    return csv_data;
}
