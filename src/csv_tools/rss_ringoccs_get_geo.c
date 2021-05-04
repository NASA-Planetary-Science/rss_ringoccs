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

rssringoccs_GeoCSV *rssringoccs_Get_Geo(const char *filename,
                                        rssringoccs_Bool use_deprecated)
{
    rssringoccs_GeoCSV *geo;
    FILE *fp;
    char buffer[1024];
    char *record, *line;
    int ch;
    unsigned long line_count, column_count, n;

    geo = malloc(sizeof(*geo));

    /*  Check if malloc failed.                                               */
    if (geo == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Get_Geo\n\n"
             "Malloc failed and returned NULL for geo. Returning.\n");
        return NULL;
    }

    /*  Initialize the pointers in the geo struct to NULL. The function       *
     *  rssringoccs_Destroy_GeoCSV_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    geo->t_oet_spm_vals = NULL;
    geo->t_ret_spm_vals = NULL;
    geo->t_set_spm_vals = NULL;
    geo->rho_km_vals = NULL;
    geo->phi_rl_deg_vals = NULL;
    geo->phi_ora_deg_vals = NULL;
    geo->B_deg_vals = NULL;
    geo->D_km_vals = NULL;
    geo->rho_dot_kms_vals = NULL;
    geo->phi_rl_dot_kms_vals = NULL;
    geo->F_km_vals = NULL;
    geo->R_imp_km_vals = NULL;
    geo->rx_km_vals = NULL;
    geo->ry_km_vals = NULL;
    geo->rz_km_vals = NULL;
    geo->vx_kms_vals = NULL;
    geo->vy_kms_vals = NULL;
    geo->vz_kms_vals = NULL;
    geo->obs_spacecract_lat_deg_vals = NULL;
    geo->error_message = NULL;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return geo;
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
    geo->n_elements = line_count;

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
    if ((column_count != 18) && (use_deprecated))
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "use_deprecated is set to true but the input CSV does not have\n"
            "18 columns. Aborting computation.\n"
        );
        return geo;
    }

    /*  And if use_deprecated is false, we need 19 column. Check this.        */
    else if ((column_count != 19) && (!use_deprecated))
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "use_deprecated is set to false but the input CSV does not have\n"
            "19 columns. Aborting computation.\n"
        );
        return geo;
    }

    /*  Allocate memory for t_oet_spm_vals and check for error.               */
    geo->t_oet_spm_vals = malloc(sizeof(*geo->t_oet_spm_vals) * line_count);
    if (geo->t_oet_spm_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_oet_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for t_ret_spm_vals and check for error.               */
    geo->t_ret_spm_vals = malloc(sizeof(*geo->t_ret_spm_vals) * line_count);
    if (geo->t_ret_spm_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_ret_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for t_set_spm_vals and check for error.               */
    geo->t_set_spm_vals = malloc(sizeof(*geo->t_set_spm_vals) * line_count);
    if (geo->t_set_spm_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_set_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for rho_km_vals and check for error.                  */
    geo->rho_km_vals = malloc(sizeof(*geo->rho_km_vals) * line_count);
    if (geo->rho_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for phi_rl_deg_vals and check for error.              */
    geo->phi_rl_deg_vals = malloc(sizeof(*geo->phi_rl_deg_vals) * line_count);
    if (geo->phi_rl_deg_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_rl_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for phi_ora_deg_vals and check for error.             */
    geo->phi_ora_deg_vals = malloc(sizeof(*geo->phi_ora_deg_vals) * line_count);
    if (geo->phi_ora_deg_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_ora_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for B_deg_vals and check for error.                   */
    geo->B_deg_vals = malloc(sizeof(*geo->B_deg_vals) * line_count);
    if (geo->B_deg_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "B_deg_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for D_km_vals and check for error.                    */
    geo->D_km_vals = malloc(sizeof(*geo->D_km_vals) * line_count);
    if (geo->D_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "D_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for rho_dot_kms_vals and check for error.             */
    geo->rho_dot_kms_vals = malloc(sizeof(*geo->rho_dot_kms_vals) * line_count);
    if (geo->rho_dot_kms_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rho_dot_kms_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for phi_rl_dot_kms_vals and check for error.          */
    geo->phi_rl_dot_kms_vals =
        malloc(sizeof(*geo->phi_rl_dot_kms_vals)*line_count);
    if (geo->phi_rl_dot_kms_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "phi_rl_dot_kms_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for F_km_vals and check for error.                    */
    geo->F_km_vals = malloc(sizeof(*geo->F_km_vals) * line_count);
    if (geo->F_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "F_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for R_imp_km_vals and check for error.                */
    geo->R_imp_km_vals = malloc(sizeof(*geo->R_imp_km_vals) * line_count);
    if (geo->R_imp_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "R_imp_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for rx_km_vals and check for error.                   */
    geo->rx_km_vals = malloc(sizeof(*geo->rx_km_vals) * line_count);
    if (geo->rx_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rx_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for ry_km_vals and check for error.                   */
    geo->ry_km_vals = malloc(sizeof(*geo->ry_km_vals) * line_count);
    if (geo->ry_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "ry_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for rz_km_vals and check for error.                   */
    geo->rz_km_vals = malloc(sizeof(*geo->rz_km_vals) * line_count);
    if (geo->rz_km_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "rz_km_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for vx_kms_vals and check for error.                  */
    geo->vx_kms_vals = malloc(sizeof(*geo->vx_kms_vals) * line_count);
    if (geo->vx_kms_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "vx_kms_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for vy_kms_vals and check for error.                  */
    geo->vy_kms_vals = malloc(sizeof(*geo->vy_kms_vals) * line_count);
    if (geo->vy_kms_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "vy_kms_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  Allocate memory for vz_kms_vals and check for error.                  */
    geo->vz_kms_vals = malloc(sizeof(*geo->vz_kms_vals) * line_count);
    if (geo->vz_kms_vals == NULL)
    {
        geo->error_occurred = rssringoccs_True;
        geo->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "vz_kms_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
        return geo;
    }

    /*  If we're using the older deprecated format, there are 18 columns. Set *
     *  the obs_spacecract_lat_deg_vals to NULL.                              */
    if (use_deprecated)
        geo->obs_spacecract_lat_deg_vals = NULL;
    else
    {
        /*  Allocate memory for obs_spacecract_lat_deg_vals.                  */
        geo->obs_spacecract_lat_deg_vals =
            malloc(sizeof(*geo->obs_spacecract_lat_deg_vals) * line_count);
        if (geo->obs_spacecract_lat_deg_vals == NULL)
        {
            geo->error_occurred = rssringoccs_True;
            geo->error_message = rssringoccs_strdup(
                "Error Encountered: rss_ringoccs\n"
                "\trssringoccs_Get_Geo\n\n"
                "Malloc returned NULL. Failed to allocate memory for.\n"
                "obs_spacecract_lat_deg_vals.\n"
                "Aborting computation and returning.\n"
            );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_GeoCSV_Members(geo);
            return geo;
        }
    }

    line = fgets(buffer, sizeof(buffer), fp);
    n = 0;
    while(line != NULL)
    {
        record = strtok(line, ",");
        geo->t_oet_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->t_ret_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->t_set_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->rho_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->phi_rl_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->phi_ora_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->B_deg_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->D_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->rho_dot_kms_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->phi_rl_dot_kms_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->F_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->R_imp_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->rx_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->ry_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->rz_km_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->vx_kms_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->vy_kms_vals[n] = atof(record);

        record = strtok(NULL, ",");
        geo->vz_kms_vals[n] = atof(record);

        record = strtok(NULL, ",");
        if (record != NULL)
            geo->obs_spacecract_lat_deg_vals[n] = atof(record);

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    return geo;
}
