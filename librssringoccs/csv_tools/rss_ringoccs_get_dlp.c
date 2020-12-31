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
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

rssringoccs_DLPCSV *rssringoccs_Get_DLP(const char *filename,
                                        rssringoccs_Bool use_deprecated)
{
    rssringoccs_DLPCSV *dlp;
    FILE *fp;
    char buffer[1024];
    char *record, *line;
    int ch;
    unsigned long line_count, column_count, n;

    dlp = malloc(sizeof(*dlp));

    /*  Check if malloc failed.                                               */
    if (dlp == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Get_DLP\n\n"
             "Malloc failed and returned NULL for dlp. Returning.\n");
        return NULL;
    }

    /*  Initialize the pointers in the dlp struct to NULL. The function       *
     *  rssringoccs_Destroy_DLPCSV_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    dlp->rho_km_vals = NULL;
    dlp->rho_corr_pole_km_vals = NULL;
    dlp->rho_corr_timing_km_vals = NULL;
    dlp->phi_rl_deg_vals = NULL;
    dlp->phi_ora_deg_vals = NULL;
    dlp->p_norm_vals = NULL;
    dlp->raw_tau_vals = NULL;
    dlp->phase_deg_vals = NULL;
    dlp->raw_tau_threshold_vals = NULL;
    dlp->t_oet_spm_vals = NULL;
    dlp->t_ret_spm_vals = NULL;
    dlp->t_set_spm_vals = NULL;
    dlp->B_deg_vals = NULL;
    dlp->error_message = NULL;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_DLP\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return dlp;
    }

    /*  Count the number of lines in the CSV.                                 */
    line_count = 0;
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
            line_count++;
    }
    rewind(fp);
    dlp->n_elements = line_count;

    /*  And count the number of columns.                                      */
    column_count = 0;
    line = fgets(buffer,sizeof(buffer), fp);
    record = strtok(line, ",");
    while (record != NULL)
    {
        record = strtok(NULL, ",");
        column_count++;
    }
    rewind(fp);

    /*  If use_deprecated was set to true, column_count must be 18. Check.    */
    if ((column_count != 12) && (use_deprecated))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_DLP\n\n"
            "use_deprecated is set to true but the input CSV does not have\n"
            "12 columns. Aborting computation.\n"
        );
        return dlp;
    }

    /*  And if use_deprecated is false, we need 19 column. Check this.        */
    else if ((column_count != 13) && (!use_deprecated))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "use_deprecated is set to false but the input CSV does not have\n"
            "13 columns. Aborting computation.\n"
        );
        return dlp;
    }

    /*  Allocate memory for t_oet_spm_vals and check for error.               */
    dlp->t_oet_spm_vals = malloc(sizeof(*dlp->t_oet_spm_vals) * line_count);
    if (dlp->t_oet_spm_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_oet_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for t_ret_spm_vals and check for error.               */
    dlp->t_ret_spm_vals = malloc(sizeof(*dlp->t_ret_spm_vals) * line_count);
    if (dlp->t_ret_spm_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_ret_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for t_set_spm_vals and check for error.               */
    dlp->t_set_spm_vals = malloc(sizeof(*dlp->t_set_spm_vals) * line_count);
    if (dlp->t_set_spm_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_set_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for rho_km_vals and check for error.                  */
    dlp->rho_km_vals = malloc(sizeof(*dlp->rho_km_vals) * line_count);
    if (dlp->rho_km_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for rho_corr_pole_km_vals and check for error.        */
    dlp->rho_corr_pole_km_vals
        = malloc(sizeof(*dlp->rho_corr_pole_km_vals) * line_count);
    if (dlp->rho_corr_pole_km_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_corr_pole_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for rho_corr_timing_km_vals and check for error.      */
    dlp->rho_corr_timing_km_vals
        = malloc(sizeof(*dlp->rho_corr_timing_km_vals) * line_count);
    if (dlp->rho_corr_timing_km_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_corr_timing_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for phi_rl_deg_vals and check for error.              */
    dlp->phi_rl_deg_vals = malloc(sizeof(*dlp->phi_rl_deg_vals) * line_count);
    if (dlp->phi_rl_deg_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_rl_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for phi_ora_deg_vals and check for error.             */
    dlp->phi_ora_deg_vals = malloc(sizeof(*dlp->phi_ora_deg_vals) * line_count);
    if (dlp->phi_ora_deg_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_ora_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for B_deg_vals and check for error.                   */
    dlp->B_deg_vals = malloc(sizeof(*dlp->B_deg_vals) * line_count);
    if (dlp->B_deg_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "B_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for raw_tau_vals and check for error.                 */
    dlp->raw_tau_vals = malloc(sizeof(*dlp->raw_tau_vals) * line_count);
    if (dlp->raw_tau_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "raw_tau_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for phase_deg_vals and check for error.               */
    dlp->phase_deg_vals = malloc(sizeof(*dlp->phase_deg_vals) * line_count);
    if (dlp->phase_deg_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phase_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  Allocate memory for raw_tau_threshold_vals and check for error.       */
    dlp->raw_tau_threshold_vals
        = malloc(sizeof(*dlp->raw_tau_threshold_vals) * line_count);
    if (dlp->raw_tau_threshold_vals == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "raw_tau_threshold_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_DLPCSV_Members(dlp);
        return dlp;
    }

    /*  If we're using the older deprecated format, there are 18 columns. Set *
     *  the p_norm_vals to NULL.                                              */
    if (use_deprecated)
        dlp->p_norm_vals = NULL;
    else
    {
        /*  Allocate memory for p_norm_vals.                                  */
        dlp->p_norm_vals =
            malloc(sizeof(*dlp->p_norm_vals) * line_count);
        if (dlp->p_norm_vals == NULL)
        {
            dlp->error_occurred = rssringoccs_True;
            dlp->error_message = rssringoccs_strdup(
                "Error Encountered: rss_ringoccs\n"
                "\trssringoccs_Get_DLP\n\n"
                "Malloc returned NULL. Failed to allocate memory for.\n"
                "p_norm_vals. Aborting computation and returning.\n"
            );

            /*  Free the variables that have been malloc'd so far.            */
            rssringoccs_Destroy_DLPCSV_Members(dlp);
            return dlp;
        }
    }

    line = fgets(buffer, sizeof(buffer), fp);
    n = 0;
    while(line != NULL)
    {
        record = strtok(line, ",");
        dlp->rho_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->rho_corr_pole_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->rho_corr_timing_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->phi_rl_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->phi_ora_deg_vals[n] = atof(record);

        if (!use_deprecated)
        {
            record = strtok(NULL, ",");
            dlp->p_norm_vals[n] = atof(record);
        }

        record = strtok(NULL, ",");
        dlp->raw_tau_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->phase_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->raw_tau_threshold_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->t_oet_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->t_ret_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->t_set_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        dlp->B_deg_vals[n] = atof(record);

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    return dlp;
}
