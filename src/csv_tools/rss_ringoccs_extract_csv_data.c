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
 *      Extracts the data from CSV files to mimic a DLP object. This object   *
 *      can then be used for diffraction reconstruction.                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Booleans, interpolation, math routines, and more.                         */
#include <libtmpl/include/tmpl.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  free and malloc are found here.                                           */
#include <stdlib.h>

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

/*  Extracts all data from CSV (.TAB) files and creates a DLP-like struct.    */
rssringoccs_CSVData *
rssringoccs_Extract_CSV_Data(const char *geo,
                             const char *cal,
                             const char *dlp,
                             const char *tau,
                             tmpl_Bool use_deprecated)
{
    /*  Pointers for the Geo, Cal, DLP, and Tau CSV objects.                  */
    rssringoccs_GeoCSV *geo_dat;
    rssringoccs_DLPCSV *dlp_dat;
    rssringoccs_CalCSV *cal_dat;
    rssringoccs_TauCSV *tau_dat;

    /*  Buffer for error messages.                                            */
    char err_mes[512];

    /*  A pointer to the CSV object.                                          */
    rssringoccs_CSVData *csv_data;

    /*  Variable for indexing.                                                */
    size_t n;

    /*  Constant for zero of type "size_t".                                   */
    const size_t zero = (size_t)0;

    /*  Variables for checking the geometry of the occultation.               */
    double min_dr_dt, max_dr_dt, temp, mu, log_power;
    double *cal_f_sky_hz_vals;

    /*  Allocate memory for the CSV object.                                   */
    csv_data = malloc(sizeof(*csv_data));

    /*  Check if malloc failed.                                               */
    if (csv_data == NULL)
        return NULL;

    /*  Initialize the members to NULL. This will prevent functions from      *
     *  trying to free pointers that weren't malloc'd in the event of error.  */
    csv_data->B_deg_vals = NULL;
    csv_data->D_km_vals = NULL;
    csv_data->f_sky_hz_vals = NULL;
    csv_data->p_norm_vals = NULL;
    csv_data->raw_tau_vals = NULL;
    csv_data->phase_deg_vals = NULL;
    csv_data->phi_deg_vals = NULL;
    csv_data->phi_rl_deg_vals = NULL;
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
    csv_data->tau_phase = NULL;
    csv_data->tau_power = NULL;
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
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned NULL for geo_dat. Aborting.\n"
        );
        return csv_data;
    }

    /*  If the CSV is empty there is something wrong with the geo string.     */
    else if (geo_dat->n_elements == zero)
    {
        csv_data->error_occurred = tmpl_True;

        if (geo_dat->error_message)
        {
            sprintf(err_mes,
                    "\nError Encountered: rss_ringoccs\n"
                    "\trssringoccs_Extract_CSV_Data\n\n"
                    "rssringoccs_Get_Geo returned an empty struct. Aborting.\n"
                    "rssringoccs_Get_Geo set the following message:\n\n%s",
                    geo_dat->error_message);
            csv_data->error_message = tmpl_strdup(err_mes);
        }
        else
        {
            csv_data->error_message = tmpl_strdup(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_Extract_CSV_Data\n\n"
                "rssringoccs_Get_Geo returned an empty struct. Aborting.\n"
            );
        }

        /*  The Geo object was allocated memory. Free before aborting.        */
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
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned NULL for dlp_dat. Aborting.\n"
        );

        /*  The Geo object was successfully created. Destroy it to avoid      *
         *  memory leaks and then abort the computation.                      */
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }

    /*  Similar check as geo, make sure the DLP object isn't empty.           */
    else if (dlp_dat->n_elements == zero)
    {
        csv_data->error_occurred = tmpl_True;

        if (dlp_dat->error_message)
        {
            sprintf(err_mes,
                    "\nError Encountered: rss_ringoccs\n"
                    "\trssringoccs_Extract_CSV_Data\n\n"
                    "rssringoccs_Get_DLP returned an empty struct. Aborting.\n"
                    "rssringoccs_Get_DLP set the following message:\n\n%s",
                    dlp_dat->error_message);
            csv_data->error_message = tmpl_strdup(err_mes);
        }
        else
        {
            csv_data->error_message = tmpl_strdup(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_Extract_CSV_Data\n\n"
                "rssringoccs_Get_DLP returned an empty struct. Aborting.\n"
            );
        }

        /*  Both Geo and DLP have memory allocated to them. Free them.        */
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
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned NULL for cal_dat. Aborting.\n"
        );

        /*  Geo and DLP were successfully created. They need to be freed to   *
         *  avoid memory leaks. Destroy both objects and then abort.          */
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }

    /*  Ensure the created object isn't empty.                                */
    else if (cal_dat->n_elements == zero)
    {
        csv_data->error_occurred = tmpl_True;

        if (cal_dat->error_message)
        {
            sprintf(err_mes,
                    "\nError Encountered: rss_ringoccs\n"
                    "\trssringoccs_Extract_CSV_Data\n\n"
                    "rssringoccs_Get_Cal returned an empty struct. Aborting.\n"
                    "rssringoccs_Get_Cal set the following message:\n\n%s",
                    cal_dat->error_message);
            csv_data->error_message = tmpl_strdup(err_mes);
        }
        else
        {
            csv_data->error_message = tmpl_strdup(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_Extract_CSV_Data\n\n"
                "rssringoccs_Get_Cal returned an empty struct. Aborting.\n"
            );
        }

        /*  Free the three CSV objects before aborting.                       */
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }

    /*  Steal the pointer to avoid a call to malloc. It is destroyed in the   *
     *  end anyways, so no harm done.                                         */
    cal_f_sky_hz_vals = cal_dat->f_sky_pred_vals;

    /*  The sky frequency is the difference of the predicted and residual     *
     *  frequencies. Loop through the arrays and compute the difference.      */
    for (n = zero; n < cal_dat->n_elements; ++n)
        cal_f_sky_hz_vals[n] -= cal_dat->f_sky_resid_fit_vals[n];

    /*  Extract the data from the TAU.TAB file.                               */
    if (tau)
    {
        tau_dat = rssringoccs_Get_Tau(tau, use_deprecated);

        /*  Check for errors.                                                 */
        if (tau_dat == NULL)
        {
            csv_data->error_occurred = tmpl_True;
            csv_data->error_message = tmpl_strdup(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_Extract_CSV_Data\n\n"
                "rssringoccs_Get_Tau returned NULL for tau_dat. Aborting.\n"
            );

            /*  Free all data before aborting.                                */
            rssringoccs_Destroy_GeoCSV(&geo_dat);
            rssringoccs_Destroy_DLPCSV(&dlp_dat);
            rssringoccs_Destroy_CalCSV(&cal_dat);
            return csv_data;
        }

        /*  Ensure the Tau object isn't empty.                                */
        else if (tau_dat->n_elements == zero)
        {
            csv_data->error_occurred = tmpl_True;

            if (tau_dat->error_message)
            {
                sprintf(
                    err_mes,
                    "\nError Encountered: rss_ringoccs\n"
                    "\trssringoccs_Extract_CSV_Data\n\n"
                    "rssringoccs_Get_Tau returned an empty struct. Aborting.\n"
                    "rssringoccs_Get_Tau set the following message:\n\n%s",
                    tau_dat->error_message
                );
                csv_data->error_message = tmpl_strdup(err_mes);
            }
            else
            {
                csv_data->error_message = tmpl_strdup(
                    "\nError Encountered: rss_ringoccs\n"
                    "\trssringoccs_Extract_CSV_Data\n\n"
                    "rssringoccs_Get_Tau returned an empty struct. Aborting.\n"
                );
            }

            /*  Geo, Cal, and DLP were successfully created, and Tau has      *
             *  memory allocated. Free all memory to avoid leaks.             */
            rssringoccs_Destroy_GeoCSV(&geo_dat);
            rssringoccs_Destroy_DLPCSV(&dlp_dat);
            rssringoccs_Destroy_CalCSV(&cal_dat);
            rssringoccs_Destroy_TauCSV(&tau_dat);
            return csv_data;
        }
    }

    /*  If Tau data is not to be extracted, avoid free'ing non-malloced       *
     *  memory by setting this to NULL. The "destroy" function will not       *
     *  attempt to free a NULL pointer.                                       */
    else
        tau_dat = NULL;

    /*  Grab the number of elements from the DLP CSV. This will be the number *
     *  of elements in the output.                                            */
    csv_data->n_elements = dlp_dat->n_elements;

    /*  The MALLOC_CSV_VAR macro contains an if-then statement with           *
     *  braces {} hence there is no need for a semi-colon at the end.         */
    MALLOC_CSV_VAR(D_km_vals)
    MALLOC_CSV_VAR(f_sky_hz_vals)
    MALLOC_CSV_VAR(rho_dot_kms_vals)
    MALLOC_CSV_VAR(rx_km_vals)
    MALLOC_CSV_VAR(ry_km_vals)
    MALLOC_CSV_VAR(rz_km_vals)

    /*  If Tau data is to be extracted, reserve memory for the variables.     */
    if (tau)
    {
        MALLOC_CSV_VAR(tau_phase)
        MALLOC_CSV_VAR(tau_power)
        MALLOC_CSV_VAR(tau_vals)
    }

    /*  The references to these variables can be stolen to avoid copying.     */
    csv_data->rho_km_vals = dlp_dat->rho_km_vals;
    csv_data->raw_tau_vals = dlp_dat->raw_tau_vals;
    csv_data->t_ret_spm_vals = dlp_dat->t_ret_spm_vals;
    csv_data->t_set_spm_vals = dlp_dat->t_set_spm_vals;
    csv_data->t_oet_spm_vals = dlp_dat->t_oet_spm_vals;
    csv_data->rho_corr_pole_km_vals = dlp_dat->rho_corr_pole_km_vals;
    csv_data->rho_corr_timing_km_vals = dlp_dat->rho_corr_timing_km_vals;
    csv_data->raw_tau_threshold_vals = dlp_dat->raw_tau_threshold_vals;

    /*  The new TAB files have normalized power included.                     */
    if (!use_deprecated)
        csv_data->p_norm_vals = dlp_dat->p_norm_vals;

    /*  Otherwise we need to compute the power from the optical depth.        */
    else
    {
        /*  Need to allocate memory for p_norm_vals.                          */
        MALLOC_CSV_VAR(p_norm_vals)

        /*  Loop through the elements of the DLP tau variable and compute the *
         *  normalize diffracted power.                                       */
        for (n = zero; n < csv_data->n_elements; ++n)
        {
            /*  tmpl_Double_Sind computes sine of an argument in degrees.     */
            mu = tmpl_Double_Sind(tmpl_Double_Abs(dlp_dat->B_deg_vals[n]));
            log_power = -csv_data->raw_tau_vals[n] / mu;

            /*  Normalize diffracted power can be computed from optical depth.*/
            csv_data->p_norm_vals[n] = tmpl_Double_Exp(log_power);
        }
    }

    /*  Various angles from the dlp file. Steal the reference.                */
    csv_data->phase_deg_vals = dlp_dat->phase_deg_vals;
    csv_data->phi_deg_vals = dlp_dat->phi_ora_deg_vals;
    csv_data->B_deg_vals = dlp_dat->B_deg_vals;
    csv_data->phi_rl_deg_vals = dlp_dat->phi_rl_deg_vals;

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
        tmpl_Double_Array_Reverse(geo_dat->rho_km_vals, geo_dat->n_elements);
        tmpl_Double_Array_Reverse(geo_dat->rho_dot_kms_vals,
                                  geo_dat->n_elements);
        tmpl_Double_Array_Reverse(geo_dat->D_km_vals, geo_dat->n_elements);
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

    tmpl_Double_Sorted_Interp1d(cal_dat->t_oet_spm_vals,
                                cal_f_sky_hz_vals,
                                cal_dat->n_elements,
                                csv_data->t_oet_spm_vals,
                                csv_data->f_sky_hz_vals,
                                csv_data->n_elements);

    /*  Interpolate the Tau data if requested.                                */
    if (tau)
    {
        if (use_deprecated)
        {
            tmpl_Double_Sorted_Interp1d(tau_dat->rho_km_vals,
                                        tau_dat->tau_vals,
                                        tau_dat->n_elements,
                                        csv_data->rho_km_vals,
                                        csv_data->tau_power,
                                        csv_data->n_elements);

            for (n = zero; n < csv_data->n_elements; ++n)
            {
                mu = tmpl_Double_Sind(tmpl_Double_Abs(csv_data->B_deg_vals[n]));
                log_power = -csv_data->tau_power[n] / mu;
                csv_data->tau_power[n] = tmpl_Double_Exp(log_power);
            }
        }
        else
            tmpl_Double_Sorted_Interp1d(tau_dat->rho_km_vals,
                                        tau_dat->power_vals,
                                        tau_dat->n_elements,
                                        csv_data->rho_km_vals,
                                        csv_data->tau_power,
                                        csv_data->n_elements);


        tmpl_Double_Sorted_Interp1d(tau_dat->rho_km_vals,
                                    tau_dat->phase_deg_vals,
                                    tau_dat->n_elements,
                                    csv_data->rho_km_vals,
                                    csv_data->tau_phase,
                                    csv_data->n_elements);

        tmpl_Double_Sorted_Interp1d(tau_dat->rho_km_vals,
                                    tau_dat->tau_vals,
                                    tau_dat->n_elements,
                                    csv_data->rho_km_vals,
                                    csv_data->tau_vals,
                                    csv_data->n_elements);
    }

    /*  Free the Geo, Cal, and TAU structs. The CSV struct stole several      *
     *  pointers from the DLP struct, so do not destroy this.                 */
    rssringoccs_Destroy_GeoCSV(&geo_dat);
    rssringoccs_Destroy_CalCSV(&cal_dat);
    rssringoccs_Destroy_TauCSV(&tau_dat);

    /*  The members of the DLP struct cannot be freed, the CSV struct has the *
     *  same pointers. We can, however, free the pointer itself. We will not  *
     *  lose references to its members since csv_data has these.              */
    free(dlp_dat);
    return csv_data;
}
