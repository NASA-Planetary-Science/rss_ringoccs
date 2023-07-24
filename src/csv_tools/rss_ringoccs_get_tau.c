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

/*  Booleans, interpolation, math routines, and more.                         */
#include <libtmpl/include/tmpl.h>

/*  Typedefs for CSV structs and function prototype given here.               */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  malloc and free are found here.                                           */
#include <stdlib.h>

/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  strtok function provided here.                                            */
#include <string.h>

rssringoccs_TauCSV *rssringoccs_Get_Tau(const char *filename,
                                        tmpl_Bool use_deprecated)
{
    rssringoccs_TauCSV *tau;
    FILE *fp;
    char buffer[1024];
    char *record, *line;
    size_t column_count, n;

    tau = malloc(sizeof(*tau));

    /*  Check if malloc failed.                                               */
    if (tau == NULL)
        return NULL;

    /*  Initialize the pointers in the tau struct to NULL. The function       *
     *  rssringoccs_Destroy_TauCSV_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    tau->rho_km_vals = NULL;
    tau->rho_corr_pole_km_vals = NULL;
    tau->rho_corr_timing_km_vals = NULL;
    tau->phi_rl_deg_vals = NULL;
    tau->phi_ora_deg_vals = NULL;
    tau->power_vals = NULL;
    tau->tau_vals = NULL;
    tau->phase_deg_vals = NULL;
    tau->tau_threshold_vals = NULL;
    tau->t_oet_spm_vals = NULL;
    tau->t_ret_spm_vals = NULL;
    tau->t_set_spm_vals = NULL;
    tau->B_deg_vals = NULL;
    tau->error_message = NULL;
    tau->error_occurred = tmpl_False;
    tau->n_elements = (size_t)0;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Tau\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return tau;
    }

    /*  And count the number of columns.                                      */
    column_count = 0;
    line = fgets(buffer,sizeof(buffer), fp);
    record = strtok(line, ",");
    while (record != NULL)
    {
        record = strtok(NULL, ",");
        column_count++;
    }

    /*  If use_deprecated was set to true, column_count must be 12. Check.    */
    if ((column_count != 12) && (use_deprecated))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Tau\n\n"
            "use_deprecated is set to true but the input CSV does not have\n"
            "12 columns. Aborting computation.\n"
        );
        fclose(fp);
        return tau;
    }

    /*  And if use_deprecated is false, we need 13 column. Check this.        */
    else if ((column_count != 13) && (!use_deprecated))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "use_deprecated is set to false but the input CSV does not have\n"
            "13 columns. Aborting computation.\n"
        );
        fclose(fp);
        return tau;
    }

    /*  Count the number of lines in the CSV.                                 */
    rewind(fp);

    tau->n_elements = tmpl_Line_Count(fp);

    /*  Allocate memory for t_oet_spm_vals and check for error.               */
    tau->t_oet_spm_vals = malloc(sizeof(*tau->t_oet_spm_vals)*tau->n_elements);
    if (tau->t_oet_spm_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_oet_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for t_ret_spm_vals and check for error.               */
    tau->t_ret_spm_vals = malloc(sizeof(*tau->t_ret_spm_vals)*tau->n_elements);
    if (tau->t_ret_spm_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_ret_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for t_set_spm_vals and check for error.               */
    tau->t_set_spm_vals = malloc(sizeof(*tau->t_set_spm_vals)*tau->n_elements);
    if (tau->t_set_spm_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_set_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for rho_km_vals and check for error.                  */
    tau->rho_km_vals = malloc(sizeof(*tau->rho_km_vals)*tau->n_elements);
    if (tau->rho_km_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for rho_corr_pole_km_vals and check for error.        */
    tau->rho_corr_pole_km_vals
        = malloc(sizeof(*tau->rho_corr_pole_km_vals)*tau->n_elements);
    if (tau->rho_corr_pole_km_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_corr_pole_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for rho_corr_timing_km_vals and check for error.      */
    tau->rho_corr_timing_km_vals
        = malloc(sizeof(*tau->rho_corr_timing_km_vals)*tau->n_elements);
    if (tau->rho_corr_timing_km_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_corr_timing_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for phi_rl_deg_vals and check for error.              */
    tau->phi_rl_deg_vals = malloc(sizeof(*tau->phi_rl_deg_vals)*tau->n_elements);
    if (tau->phi_rl_deg_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_rl_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for phi_ora_deg_vals and check for error.             */
    tau->phi_ora_deg_vals = malloc(sizeof(*tau->phi_ora_deg_vals) * tau->n_elements);
    if (tau->phi_ora_deg_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_ora_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for B_deg_vals and check for error.                   */
    tau->B_deg_vals = malloc(sizeof(*tau->B_deg_vals) * tau->n_elements);
    if (tau->B_deg_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "B_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for tau_vals and check for error.                     */
    tau->tau_vals = malloc(sizeof(*tau->tau_vals) * tau->n_elements);
    if (tau->tau_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "tau_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for phase_deg_vals and check for error.               */
    tau->phase_deg_vals = malloc(sizeof(*tau->phase_deg_vals) * tau->n_elements);
    if (tau->phase_deg_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phase_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  Allocate memory for tau_threshold_vals and check for error.       */
    tau->tau_threshold_vals
        = malloc(sizeof(*tau->tau_threshold_vals) * tau->n_elements);
    if (tau->tau_threshold_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Tau\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "tau_threshold_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_TauCSV_Members(tau);
        return tau;
    }

    /*  If we're using the older deprecated format, there are 18 columns. Set *
     *  the power_vals to NULL.                                               */
    if (use_deprecated)
        tau->power_vals = NULL;
    else
    {
        /*  Allocate memory for power_vals.                                   */
        tau->power_vals =
            malloc(sizeof(*tau->power_vals) * tau->n_elements);
        if (tau->power_vals == NULL)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "Error Encountered: rss_ringoccs\n"
                "\trssringoccs_Get_Tau\n\n"
                "Malloc returned NULL. Failed to allocate memory for.\n"
                "power_vals. Aborting computation and returning.\n"
            );

            /*  Free the variables that have been malloc'd so far.            */
            rssringoccs_Destroy_TauCSV_Members(tau);
            return tau;
        }
    }

    line = fgets(buffer, sizeof(buffer), fp);
    n = 0;
    while(line != NULL)
    {
        record = strtok(line, ",");
        tau->rho_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->rho_corr_pole_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->rho_corr_timing_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->phi_rl_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->phi_ora_deg_vals[n] = atof(record);

        if (!use_deprecated)
        {
            record = strtok(NULL, ",");
            tau->power_vals[n] = atof(record);
        }

        record = strtok(NULL, ",");
        tau->tau_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->phase_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->tau_threshold_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->t_oet_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->t_ret_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->t_set_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        tau->B_deg_vals[n] = atof(record);

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    return tau;
}
