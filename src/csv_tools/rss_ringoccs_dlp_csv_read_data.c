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
 *      Function for allocating memory for a DLP CSV object.                  *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
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

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for reading data from a DLP CSV into a DLP object.               */
void rssringoccs_DLPCSV_Read_Data(rssringoccs_DLPCSV *dlp, FILE *fp)
{
    /*  Buffer for reading the characters in a line from the CSV.             */
    char buffer[1024];

    /*  Variable for keeping track of the current column being read.          */
    char *record;

    /*  Pointer used for storing the current line being read in the CSV.      */
    char *line;

    /*  The constant zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  Index for the current row being read into the DLP object.             */
    size_t n = zero;

    /*  If the input DLP object is NULL there is nothing to do. Return.       */
    if (!dlp)
        return;

    /*  Similarly if an error already occurred. Abort the computation.        */
    if (dlp->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_DLPCSV_Read_Data\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (dlp->n_elements == zero)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_DLPCSV_Read_Data\n\n"
            "n_elements is zero, nothing to read. Aborting.\n"
        );

        return;
    }

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);

    while(line != NULL)
    {
        /*  The order in the CSV is:                                          *
         *      ring radius:            rho_km_vals                           *
         *      pole corrected radius:  rho_corr_pole_km_vals                 *
         *      time corrected radius:  rho_corr_timing_km_vals               *
         *      ring longitude angle:   phi_rl_deg_vals                       *
         *      ring azimuth angle:     phi_ora_deg_vals                      *
         *      normalized power:       p_norm_vals                           *
         *          Only if not using                                         *
         *          the deprecated format.                                    *
         *      optical depth:          raw_tau_vals                          *
         *      complex phase:          phase_deg_vals                        *
         *      depth threshold:        raw_tau_threshold_vals                *
         *      observed time:          t_oet_spm_vals                        *
         *      ring time:              t_ret_spm_vals                        *
         *      spacecraft time:        t_set_spm_vals                        *
         *      opening angle:          B_deg_vals                            *
         *  Read them from the CSV and cast them to double using atof.        */
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

        if (!dlp->use_deprecated)
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

    /*  Reset the file pointer back to the start of the file.                 */
    rewind(fp);
}
/*  End of rssringoccs_DLPCSV_Read_Data.                                      */
