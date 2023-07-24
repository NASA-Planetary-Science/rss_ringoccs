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
 *  Purpose:                                                                  *
 *      Extracts all geometry data from a Geo CSV.                            *
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

/*  Check if the macro name is available.                                     */
#ifdef MALLOC_GEO_VAR
#undef MALLOC_GEO_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_GEO_VAR(var)                                                    \
    geo->var = malloc(sizeof(*geo->var) * geo->n_elements);                    \
    if (geo->var == NULL)                                                      \
    {                                                                          \
        geo->error_occurred = tmpl_True;                                       \
        geo->error_message = tmpl_strdup(                                      \
            "Error Encountered: rss_ringoccs\n"                                \
            "\ttrssringoccs_Get_Cal\n\n"                                       \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_Destroy_GeoCSV_Members(geo);                               \
        fclose(fp);                                                            \
        return geo;                                                            \
    }

/*  Function for extracting the data from a CAL.TAB file.                     */
rssringoccs_GeoCSV *
rssringoccs_Get_Geo(const char *filename, tmpl_Bool use_deprecated)
{
    /*  Pointer to the Geo struct.                                            */
    rssringoccs_GeoCSV *geo;

    /*  File pointer for the CSV file.                                        */
    FILE *fp;

    /*  Buffer for reading the file line by line.                             */
    char buffer[1024];

    /*  Variables for parsing the contents of the CSV file.                   */
    char *record, *line;

    unsigned int column_count = 0U;
    size_t n = 0;

    /*  Allocate memory for the geo data.                                     */
    geo = malloc(sizeof(*geo));

    /*  Check if malloc failed.                                               */
    if (geo == NULL)
        return NULL;

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
    geo->n_elements = 0UL;
    geo->error_occurred = tmpl_False;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        geo->error_occurred = tmpl_True;
        geo->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return geo;
    }

    /*  Count the number of columns.                                          */
    line = fgets(buffer, sizeof(buffer), fp);

    /*  GEO.TAB files are comma separeted, so count the number of commas.     */
    record = strtok(line, ",");

    while (record != NULL)
    {
        record = strtok(NULL, ",");
        column_count++;

        /*  The GEO file should have either 18 or 19 columns.                 */
        if (column_count > 19U)
            break;
    }

    /*  If use_deprecated was set to true, column_count must be 18. Check.    */
    if ((column_count != 18U) && (use_deprecated))
    {
        geo->error_occurred = tmpl_True;
        geo->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "use_deprecated is set to true but the input CSV does not have\n"
            "18 columns. Aborting computation.\n"
        );
        fclose(fp);
        return geo;
    }

    /*  And if use_deprecated is false, we need 19 column. Check this.        */
    else if ((column_count != 19U) && (!use_deprecated))
    {
        geo->error_occurred = tmpl_True;
        geo->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Geo\n\n"
            "use_deprecated is set to false but the input CSV does not have\n"
            "19 columns. Aborting computation.\n"
        );
        fclose(fp);
        return geo;
    }

    /*  Reset the file back to the start.                                     */
    rewind(fp);

    /*  Count the number of lines in the CSV.                                 */
    geo->n_elements = tmpl_Line_Count(fp);

    /*  Use the MALLOC_GEO_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_GEO_VAR(t_oet_spm_vals)
    MALLOC_GEO_VAR(t_ret_spm_vals)
    MALLOC_GEO_VAR(t_set_spm_vals)
    MALLOC_GEO_VAR(rho_km_vals)
    MALLOC_GEO_VAR(phi_rl_deg_vals)
    MALLOC_GEO_VAR(phi_ora_deg_vals)
    MALLOC_GEO_VAR(B_deg_vals)
    MALLOC_GEO_VAR(D_km_vals)
    MALLOC_GEO_VAR(rho_dot_kms_vals)
    MALLOC_GEO_VAR(phi_rl_dot_kms_vals)
    MALLOC_GEO_VAR(F_km_vals)
    MALLOC_GEO_VAR(R_imp_km_vals)
    MALLOC_GEO_VAR(rx_km_vals)
    MALLOC_GEO_VAR(ry_km_vals)
    MALLOC_GEO_VAR(rz_km_vals)
    MALLOC_GEO_VAR(vx_kms_vals)
    MALLOC_GEO_VAR(vy_kms_vals)
    MALLOC_GEO_VAR(vz_kms_vals)

    /*  If we're using the older deprecated format, there are 19 columns.     *
     *  Allocate memory for obs_spacecract_lat_deg_vals.                      */
    if (!use_deprecated)
    {
        MALLOC_GEO_VAR(obs_spacecract_lat_deg_vals)
    }

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);

    while (line != NULL)
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

        if (!use_deprecated)
        {
            record = strtok(NULL, ",");
            geo->obs_spacecract_lat_deg_vals[n] = atof(record);
        }

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    /*  Close the file.                                                       */
    fclose(fp);
    return geo;
}
/*  End of rssringoccs_Get_GEO.                                               */

/*  Undefine the Macro function.                                              */
#undef MALLOC_GEO_VAR
