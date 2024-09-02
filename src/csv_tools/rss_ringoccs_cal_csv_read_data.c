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
 *      Function for reading data from a CSV file into a Cal CSV object.      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  fgets, rewind, and NULL found here.                                       */
#include <stdio.h>

/*  strtok provided here.                                                     */
#include <string.h>

/*  atof function found here.                                                 */
#include <stdlib.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  rssringoccs_CalCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for reading data from a Cal CSV into a Cal object.               */
void rssringoccs_CalCSV_Read_Data(rssringoccs_CalCSV *cal, FILE *fp)
{
    /*  Buffer for reading the characters in a line from the CSV.             */
    char buffer[1024];

    /*  Variable for keeping track of the current column being read.          */
    char *record;

    /*  Pointer used for storing the current line being read in the CSV.      */
    char *line;

    /*  The constant zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  Index for the current row being read into the Cal object.             */
    size_t n = zero;

    /*  If the input Cal object is NULL there is nothing to do. Return.       */
    if (!cal)
        return;

    /*  Similarly if an error already occurred. Abort the computation.        */
    if (cal->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Read_Data\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (cal->n_elements == zero)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Read_Data\n\n"
            "n_elements is zero, nothing to read. Aborting.\n"
        );

        return;
    }

    /*  Read in the data and start copying it to the Cal object.              */
    line = fgets(buffer, sizeof(buffer), fp);

    while (line != NULL)
    {
        /*  The order in the CSV is:                                          *
         *      time:                   t_oet_spm_vals                        *
         *      frequency (predicted):  f_sky_pred_vals                       *
         *      frequency (residual):   f_sky_resid_fit_vals                  *
         *      power (free-space):     p_free_vals                           *
         *  Read them from the CSV and cast them to double using atof.        */
        record = strtok(line, ",");
        cal->t_oet_spm_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->f_sky_pred_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->f_sky_resid_fit_vals[n] = atof(record);

        record = strtok(NULL, ",");
        cal->p_free_vals[n] = atof(record);

        /*  Read the next row in the CSV and increment the index.             */
        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    /*  Reset the file pointer back to the start of the file.                 */
    rewind(fp);
}
/*  End of rssringoccs_CalCSV_Read_Data.                                      */
