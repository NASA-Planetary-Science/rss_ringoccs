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

#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

rssringoccs_CalCSV *rssringoccs_Get_Cal(const char *filename)
{
    rssringoccs_CalCSV *cal;
    FILE *fp;
    char buffer[1024];
    char *record, *line;
    int ch;
    unsigned long line_count, column_count, n;

    cal = malloc(sizeof(*cal));

    /*  Check if malloc failed.                                               */
    if (cal == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Get_Cal\n\n"
             "Malloc failed and returned NULL for cal. Returning.\n");
        return NULL;
    }

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
    line_count = 0;
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
            line_count++;
    }
    rewind(fp);
    cal->n_elements = line_count;

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

    /*  Allocate memory for t_oet_spm_vals and check for error.               */
    cal->t_oet_spm_vals = malloc(sizeof(*cal->t_oet_spm_vals) * line_count);
    if (cal->t_oet_spm_vals == NULL)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "t_oet_spm_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_CalCSV_Members(cal);
        return cal;
    }

    /*  Allocate memory for f_sky_pred_vals and check for error.              */
    cal->f_sky_pred_vals = malloc(sizeof(*cal->f_sky_pred_vals) * line_count);
    if (cal->f_sky_pred_vals == NULL)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "f_sky_pred_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_CalCSV_Members(cal);
        return cal;
    }

    /*  Allocate memory for f_sky_resid_fit_vals and check for error.         */
    cal->f_sky_resid_fit_vals
        = malloc(sizeof(*cal->f_sky_resid_fit_vals) * line_count);
    if (cal->f_sky_resid_fit_vals == NULL)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "f_sky_resid_fit_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_CalCSV_Members(cal);
        return cal;
    }

    /*  Allocate memory for p_free_vals and check for error.                  */
    cal->p_free_vals = malloc(sizeof(*cal->p_free_vals) * line_count);
    if (cal->p_free_vals == NULL)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\ttrssringoccs_Get_Cal\n\n"
            "Malloc returned NULL. Failed to allocate memory for.\n"
            "p_free_vals. Aborting computation and returning.\n"
        );

        /*  Free the variables that have been malloc'd so far.                */
        rssringoccs_Destroy_CalCSV_Members(cal);
        return cal;
    }

    line = fgets(buffer, sizeof(buffer), fp);
    n = 0;
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

    return cal;
}
