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

#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*  Check if the macro name is available.                                     */
#ifdef MALLOC_GEO_VAR
#undef MALLOC_GEO_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_GEO_VAR(var)                                                    \
    cal->var = malloc(sizeof(*cal->var) * line_count);                         \
    if (cal->var == NULL)                                                      \
    {                                                                          \
        cal->error_occurred = tmpl_True;                                       \
        cal->error_message = tmpl_strdup(                                      \
            "Error Encountered: rss_ringoccs\n"                                \
            "\ttrssringoccs_Get_Cal\n\n"                                       \
            "Malloc returned NULL. Failed to allocate memory for var.\n"       \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_Destroy_CalCSV_Members(cal);                               \
        return cal;                                                            \
    }                                                                          \

/*  Function for extracting the data from a CAL.TAB file.                     */
rssringoccs_CalCSV *rssringoccs_Get_Cal(const char *filename)
{
    /*  Pointer to a CalCSV struct.                                           */
    rssringoccs_CalCSV *cal;

    /*  File object for the file we're reading.                               */
    FILE *fp;
    char buffer[1024];
    char *record, *line;
    int ch;
    unsigned int column_count;
    unsigned long int line_count, n;

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

    /*  Count the number of lines in the CSV.                                 */
    line_count = 0UL;
    while (!feof(fp))
    {
        ch = fgetc(fp);
        if (ch == '\n')
            line_count++;
    }

    /*  Reset the file back to the start.                                     */
    rewind(fp);
    cal->n_elements = line_count;

    /*  And count the number of columns.                                      */
    column_count = 0U;
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

    /*  Reset the file back to the start.                                     */
    rewind(fp);

    /*  There should be 4 columns. Check this.                                */
    if (column_count != 4)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Input CSV does not have 4 columns. Aborting computation.\n"
        );
        return cal;
    }

    /*  Use the MALLOC_GEO_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_GEO_VAR(t_oet_spm_vals)
    MALLOC_GEO_VAR(f_sky_pred_vals)
    MALLOC_GEO_VAR(f_sky_resid_fit_vals)
    MALLOC_GEO_VAR(p_free_vals)

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);
    n = 0U;
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
#undef MALLOC_GEO_VAR

