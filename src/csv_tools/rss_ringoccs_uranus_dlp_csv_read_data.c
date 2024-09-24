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
#define DLP_READ_VAR(var)                                                      \
    record = strtok(NULL, ",");                                                \
    dlp->var[n] = atof(record)

#define RSS_RINGOCCS_DEG2RAD (+5.7295779513082320876798154814105170E+01)

/*  Function for reading data from a DLP CSV into a DLP object.               */
void rssringoccs_UranusDLPCSV_Read_Data(rssringoccs_UranusDLPCSV *dlp, FILE *fp)
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
            "\trssringoccs_UranusDLPCSV_Read_Data\n\n"
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
            "\trssringoccs_UranusDLPCSV_Read_Data\n\n"
            "n_elements is zero, nothing to read. Aborting.\n"
        );

        return;
    }

    /*  Read in all of the data.                                              */
    line = fgets(buffer, sizeof(buffer), fp);

    while (line != NULL)
    {
        /*  The order in the CSV is:                                          *
         *      rho_km_vals:            Ring radius.                          *
         *      phi_rl_deg_vals:        Ring longitude angle.                 *
         *      phi_deg_vals:           Ring azimuth angle.                   *
         *      p_norm_vals:            Raw diffracted power.                 *
         *      raw_tau_vals:           Raw optical depth.                    *
         *      phase_deg_vals:         Diffracted optical phase.             *
         *      raw_tau_threshold_vals: Raw optical threshold.                *
         *      t_oet_spm_vals:         Observed event time.                  *
         *      t_ret_spm_vals:         Ring event time.                      *
         *      t_set_spm_vals:         Spacecraft event time.                *
         *      B_deg_vals:             Ring opening angle.                   *
         *      rho_dot_kms_vals        Radial velocity in ring plane.        *
         *      F_km_vals               Fresnel scale.                        *
         *      D_km_vals               Spacecraft to ring plane distance.    *
         *      f_sky_hz_vals:          Sky frequency, in Hertz.              *
         *  Read them from the CSV and cast them to double using atof.        */
        record = strtok(line, ",");
        dlp->rho_km_vals[n] = atof(record);

        /*  The DLP_READ_VAR macro looks at the next column in the CSV and    *
         *  converts the data into a double, storing it in the struct.        */
        DLP_READ_VAR(rho_corr_pole_km_vals);
        DLP_READ_VAR(rho_corr_timing_km_vals);
        DLP_READ_VAR(phi_rl_deg_vals);
        DLP_READ_VAR(phi_deg_vals);
        DLP_READ_VAR(p_norm_vals);
        DLP_READ_VAR(raw_tau_vals);
        DLP_READ_VAR(phase_deg_vals);
        DLP_READ_VAR(raw_tau_threshold_vals);
        DLP_READ_VAR(t_oet_spm_vals);
        DLP_READ_VAR(t_ret_spm_vals);
        DLP_READ_VAR(t_set_spm_vals);
        DLP_READ_VAR(B_deg_vals);
        DLP_READ_VAR(rho_dot_kms_vals);
        DLP_READ_VAR(F_km_vals);
        DLP_READ_VAR(D_km_vals);
        DLP_READ_VAR(f_sky_hz_vals);

        /*  Read the next line in the CSV and increment the index.            */
        line = fgets(buffer, sizeof(buffer), fp);
        ++n;
    }

    /*  Reset the file pointer back to the start of the file.                 */
    rewind(fp);

    /*  If the data was in radians, convert it to degrees.                    */
    if (dlp->in_radians)
    {
        for (n = 0; n < dlp->n_elements; ++n)
        {
            dlp->phi_rl_deg_vals[n] *= RSS_RINGOCCS_DEG2RAD;
            dlp->phi_deg_vals[n] *= RSS_RINGOCCS_DEG2RAD;
            dlp->phase_deg_vals[n] *= RSS_RINGOCCS_DEG2RAD;
            dlp->B_deg_vals[n] *= RSS_RINGOCCS_DEG2RAD;
        }

        /*  The data is no longer in radians. Flip this variable.             */
        dlp->in_radians = tmpl_False;
    }
}
/*  End of rssringoccs_UranusDLPCSV_Read_Data.                                */

/*  Undefine everything in case someone wants to #include this file.          */
#undef DLP_READ_VAR
#undef RSS_RINGOCCS_DEG2RAD
