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
 *      Extracts all data from a Calibration CSV.                             *
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
#ifdef MALLOC_CAL_VAR
#undef MALLOC_CAL_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_CAL_VAR(var)                                                    \
    cal->var = malloc(sizeof(*cal->var) * cal->n_elements);                    \
    if (cal->var == NULL)                                                      \
    {                                                                          \
        cal->error_occurred = tmpl_True;                                       \
        cal->error_message = tmpl_strdup(                                      \
            "Error Encountered: rss_ringoccs\n"                                \
            "\ttrssringoccs_Get_Cal\n\n"                                       \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_Destroy_CalCSV_Members(cal);                               \
        fclose(fp);                                                            \
        return cal;                                                            \
    }

/*  Function for extracting the data from a CAL.TAB file.                     */
rssringoccs_CalCSV *rssringoccs_Get_Cal(const char *filename)
{
    /*  Pointer to a CalCSV struct.                                           */
    rssringoccs_CalCSV *cal;

    /*  File object for the file we're reading.                               */
    FILE *fp;

    /*  Buffer for reading a line of the CSV into.                            */
    char buffer[1024];

    /*  Variables used for parsing the contents of the CSV.                   */
    char *record, *line;

    unsigned int column_count = 0U;
    size_t n = 0;

    /*  Allocate memory for the CalCSV object.                                */
    cal = malloc(sizeof(*cal));

    /*  Check if malloc failed.                                               */
    if (cal == NULL)
        return NULL;

    /*  Initialize the pointers in the cal struct to NULL. The function       *
     *  rssringoccs_Destroy_CalCSV_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    cal->t_oet_spm_vals = NULL;
    cal->f_sky_pred_vals = NULL;
    cal->f_sky_resid_fit_vals = NULL;
    cal->p_free_vals = NULL;
    cal->error_message = NULL;
    cal->n_elements = 0UL;
    cal->error_occurred = tmpl_False;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (fp == NULL)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Cal\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );
        return cal;
    }

    /*  Count the number of columns.                                          */
    line = fgets(buffer, sizeof(buffer), fp);

    /*  CAL.TAB files are comma separeted, so count the number of commas.     */
    record = strtok(line, ",");

    while (record != NULL)
    {
        record = strtok(NULL, ",");
        column_count++;

        /*  The CAL file should have exactly 4 columns.                       */
        if (column_count > 4U)
            break;
    }

    /*  There should be 4 columns. Check this.                                */
    if (column_count != 4U)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Input CSV does not have 4 columns. Aborting computation.\n"
        );
        fclose(fp);
        return cal;
    }

    /*  Reset the file back to the start.                                     */
    rewind(fp);

    /*  Count the number of lines in the CSV.                                 */
    cal->n_elements = tmpl_Line_Count(fp);

    /*  Use the MALLOC_CAL_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_CAL_VAR(t_oet_spm_vals)
    MALLOC_CAL_VAR(f_sky_pred_vals)
    MALLOC_CAL_VAR(f_sky_resid_fit_vals)
    MALLOC_CAL_VAR(p_free_vals)

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);

    while(line != NULL)
    {
        record = strtok(line, ",");
        cal->t_oet_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->f_sky_pred_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->f_sky_resid_fit_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->p_free_vals[n] = atof(record);

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    /*  Close the file.                                                       */
    fclose(fp);
    return cal;
}
/*  End of rssringoccs_Get_Cal.                                               */

/*  Undefine the Macro function.                                              */
#undef MALLOC_CAL_VAR
