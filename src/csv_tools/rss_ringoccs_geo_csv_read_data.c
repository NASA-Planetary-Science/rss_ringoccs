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
 *      Function for reading data from a CSV file into a Geo CSV object.      *
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

/*  libtmpl provides Booleans.                                                */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_GeoCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for reading data from a Geo CSV into a Geo object.               */
void rssringoccs_GeoCSV_Read_Data(rssringoccs_GeoCSV *geo, FILE *fp)
{
    /*  Buffer for reading the characters in a line from the CSV.             */
    char buffer[1024];

    /*  Variable for keeping track of the current column being read.          */
    char *record;

    /*  Pointer used for storing the current line being read in the CSV.      */
    char *line;

    /*  The constant zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  Index for the current row being read into the Geo object.             */
    size_t n = zero;

    /*  If the input Geo object is NULL there is nothing to do. Return.       */
    if (!geo)
        return;

    /*  Similarly if an error already occurred. Abort the computation.        */
    if (geo->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        geo->error_occurred = tmpl_True;
        geo->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_GeoCSV_Read_Data\n\n"
            "Input file is NULL.\n";

        return;
    }

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (geo->n_elements == zero)
    {
        geo->error_occurred = tmpl_True;
        geo->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_GeoCSV_Read_Data\n\n"
            "n_elements is zero, nothing to read.\n";

        return;
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

        if (!geo->use_deprecated)
        {
            record = strtok(NULL, ",");
            geo->obs_spacecraft_lat_deg_vals[n] = atof(record);
        }

        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }
    /*  Reset the file pointer back to the start of the file.                 */
    rewind(fp);
}
/*  End of rssringoccs_GeoCSV_Read_Data.                                      */
