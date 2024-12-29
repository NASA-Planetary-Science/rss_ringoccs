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
 *      Function for reading data from a CSV file into a DLP CSV object.      *
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

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Macro for reading data from the DLP CSV file and converting it to a       *
 *  double. The result is stored in the corresponding variable in the struct. */
#define DLPM_READ_VAR(var)                                                     \
    record = strtok(NULL, ",");                                                \
    dlpm->var[n] = atof(record)

/*  Function for reading data from a DLP CSV into a DLP object.               */
void
rssringoccs_MergedCSVData_Read_Data(rssringoccs_MergedCSVData *dlpm, FILE *fp)
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
    if (!dlpm)
        return;

    /*  Similarly if an error already occurred. Abort the computation.        */
    if (dlpm->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Read_Data\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (dlpm->n_elements == zero)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Read_Data\n\n"
            "n_elements is zero, nothing to read. Aborting.\n"
        );

        return;
    }

    /*  Read in the first bit of data.                                        */
    line = fgets(buffer, sizeof(buffer), fp);

    /*  Check that fgets was able to read in some data from the CSV file.     */
    if (!line)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Read_Data\n\n"
            "fgets returned NULL, no data to read. Aborting.\n"
        );

        return;
    }

    /*  If the CSV has a header, skip this. The first word of the header      *
     *  would be either 'B_deg_vals' or 'b_deg_vals'. Check for the letter B. */
    if ((line[0] == 'B') || (line[0] == 'b'))
    {
        line = fgets(buffer, sizeof(buffer), fp);

        /*  Since the starting row does not represent actual data, the number *
         *  of lines is equal to the number of elements plus one. Decrement   *
         *  the n_elements variable to set the correct value.                 */
        --dlpm->n_elements;
    }

    /*  Loop through the lines in the CSV file and read them into the struct. */
    while (line != NULL)
    {
        /*  CSV order:                                                        *
         *      B_deg_vals                                                    *
         *      D_km_vals                                                     *
         *      f_sky_hz_vals                                                 *
         *      p_norm_vals                                                   *
         *      phase_deg_vals                                                *
         *      phi_deg_vals                                                  *
         *      phi_rl_deg_vals                                               *
         *      raw_tau_threshold_vals                                        *
         *      raw_tau_vals                                                  *
         *      rho_corr_pole_km_vals                                         *
         *      rho_corr_timing_km_vals                                       *
         *      rho_dot_kms_vals                                              *
         *      rho_km_vals                                                   *
         *      rx_km_vals                                                    *
         *      ry_km_vals                                                    *
         *      rz_km_vals                                                    *
         *      t_oet_spm_vals                                                *
         *      t_ret_spm_vals                                                *
         *      t_set_spm_vals                                                *
         *  Read them from the CSV and cast them to double using atof.        */
        record = strtok(line, ",");
        dlpm->B_deg_vals[n] = atof(record);

        /*  The DLPM_READ_VAR macro looks at the next column in the CSV and   *
         *  converts the data into a double, storing it in the struct.        */
        DLPM_READ_VAR(D_km_vals);
        DLPM_READ_VAR(f_sky_hz_vals);
        DLPM_READ_VAR(p_norm_vals);
        DLPM_READ_VAR(phase_deg_vals);
        DLPM_READ_VAR(phi_deg_vals);
        DLPM_READ_VAR(phi_rl_deg_vals);
        DLPM_READ_VAR(raw_tau_threshold_vals);
        DLPM_READ_VAR(raw_tau_vals);
        DLPM_READ_VAR(rho_corr_pole_km_vals);
        DLPM_READ_VAR(rho_corr_timing_km_vals);
        DLPM_READ_VAR(rho_dot_kms_vals);
        DLPM_READ_VAR(rho_km_vals);
        DLPM_READ_VAR(rx_km_vals);
        DLPM_READ_VAR(ry_km_vals);
        DLPM_READ_VAR(rz_km_vals);
        DLPM_READ_VAR(t_oet_spm_vals);
        DLPM_READ_VAR(t_ret_spm_vals);
        DLPM_READ_VAR(t_set_spm_vals);

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    /*  Reset the file pointer back to the start of the file.                 */
    rewind(fp);
}
/*  End of rssringoccs_MergedCSVData_Read_Data.                               */

/*  Undefine everything in case someone wants to #include this file.          */
#undef DLPM_READ_VAR