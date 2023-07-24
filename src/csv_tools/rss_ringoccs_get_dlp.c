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
 *  Purpose:                                                                  *
 *      Extracts all diffraction data from a DLP CSV.                         *
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
#ifdef MALLOC_DLP_VAR
#undef MALLOC_DLP_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_DLP_VAR(var)                                                    \
    dlp->var = malloc(sizeof(*dlp->var) * dlp->n_elements);                    \
    if (dlp->var == NULL)                                                      \
    {                                                                          \
        dlp->error_occurred = tmpl_True;                                       \
        dlp->error_message = tmpl_strdup(                                      \
            "Error Encountered: rss_ringoccs\n"                                \
            "\ttrssringoccs_Get_DLP\n\n"                                       \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_Destroy_DLPCSV_Members(dlp);                               \
        fclose(fp);                                                            \
        return dlp;                                                            \
    }

/*  Function for extracting the data from a DLP.TAB file.                     */
rssringoccs_DLPCSV *
rssringoccs_Get_DLP(const char *filename, tmpl_Bool use_deprecated)
{
    /*  Pointer to the DLP struct.                                            */
    rssringoccs_DLPCSV *dlp;

    /*  File pointer for the CSV file.                                        */
    FILE *fp;

    /*  Buffer for reading the file line by line.                             */
    char buffer[1024];

    /*  Variables for parsing the contents of the CSV file.                   */
    char *record, *line;

    unsigned int column_count = 0U;
    size_t n = 0;

    /*  Allocate memory for the dlp data.                                     */
    dlp = malloc(sizeof(*dlp));

    /*  Check if malloc failed.                                               */
    if (dlp == NULL)
        return NULL;

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
    dlp->n_elements = 0UL;
    dlp->error_occurred = tmpl_False;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_DLP\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return dlp;
    }

    /*  Count the number of columns.                                          */
    line = fgets(buffer, sizeof(buffer), fp);

    /*  DLP.TAB files are comma separeted, so count the number of commas.     */
    record = strtok(line, ",");

    while (record != NULL)
    {
        record = strtok(NULL, ",");
        column_count++;

        /*  The DLP file should have 12 or 13 columns.                        */
        if (column_count > 13U)
            break;
    }

    /*  If use_deprecated was set to true, column_count must be 12. Check.    */
    if ((column_count != 12) && (use_deprecated))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_DLP\n\n"
            "use_deprecated is set to true but the input CSV does not have\n"
            "12 columns. Aborting computation.\n"
        );
        fclose(fp);
        return dlp;
    }

    /*  And if use_deprecated is false, we need 13 column. Check this.        */
    else if ((column_count != 13) && (!use_deprecated))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_DLP\n\n"
            "use_deprecated is set to false but the input CSV does not have\n"
            "13 columns. Aborting computation.\n"
        );
        fclose(fp);
        return dlp;
    }

    /*  Reset the file back to the start.                                     */
    rewind(fp);

    /*  Count the number of lines in the CSV.                                 */
    dlp->n_elements = tmpl_Line_Count(fp);

    /*  Use the MALLOC_DLP_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_DLP_VAR(t_oet_spm_vals)
    MALLOC_DLP_VAR(t_ret_spm_vals)
    MALLOC_DLP_VAR(t_set_spm_vals)
    MALLOC_DLP_VAR(rho_km_vals)
    MALLOC_DLP_VAR(rho_corr_pole_km_vals)
    MALLOC_DLP_VAR(rho_corr_timing_km_vals)
    MALLOC_DLP_VAR(phi_rl_deg_vals)
    MALLOC_DLP_VAR(phi_ora_deg_vals)
    MALLOC_DLP_VAR(B_deg_vals)
    MALLOC_DLP_VAR(raw_tau_vals)
    MALLOC_DLP_VAR(phase_deg_vals)
    MALLOC_DLP_VAR(raw_tau_threshold_vals)


    /*  If we're using the older deprecated format, there are 13 columns.     *
     *  Allocate memory for p_norm_vals as well.                              */
    if (!use_deprecated)
    {
        MALLOC_DLP_VAR(p_norm_vals)
    }

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);

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

    /*  Close the file.                                                       */
    fclose(fp);
    return dlp;
}
/*  End of rssringoccs_Get_DLP.                                               */

/*  Undefine the Macro function.                                              */
#undef MALLOC_DLP_VAR
