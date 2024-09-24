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
 *      Function for allocating memory for a Uranus DLP CSV object.           *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 24, 2024                                            *
 ******************************************************************************/

/*  malloc is found here.                                                     */
#include <stdlib.h>

/*  libtmpl provided Booleans, string duplicate, and line count.              */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if the macro name is available.                                     */
#ifdef MALLOC_DLP_VAR
#undef MALLOC_DLP_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_DLP_VAR(var)                                                    \
    dlp->var = malloc(sizeof(*dlp->var) * dlp->n_elements);                    \
    if (!dlp->var)                                                             \
    {                                                                          \
        dlp->error_occurred = tmpl_True;                                       \
        dlp->error_message = tmpl_String_Duplicate(                            \
            "Error Encountered: rss_ringoccs\n"                                \
            "\trssringoccs_UranusDLPCSV_Malloc\n\n"                            \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_UranusDLPCSV_Destroy_Members(dlp);                         \
        return;                                                                \
    }

/*  Function for allocating memory to a DLP CSV based on a CSV file pointer.  */
void rssringoccs_UranusDLPCSV_Malloc(rssringoccs_UranusDLPCSV *dlp, FILE *fp)
{
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
            "\trssringoccs_UranusDLPCSV_Malloc\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  Count the number of lines in the CSV.                                 */
    dlp->n_elements = tmpl_Line_Count(fp);

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (dlp->n_elements == (size_t)0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusDLPCSV_Malloc\n\n"
            "n_elements is zero, nothing to malloc. Aborting.\n"
        );

        return;
    }

    /*  Use the MALLOC_DLP_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_DLP_VAR(rho_km_vals)
    MALLOC_DLP_VAR(rho_corr_pole_km_vals)
    MALLOC_DLP_VAR(rho_corr_timing_km_vals)
    MALLOC_DLP_VAR(phi_rl_deg_vals)
    MALLOC_DLP_VAR(phi_deg_vals)
    MALLOC_DLP_VAR(p_norm_vals)
    MALLOC_DLP_VAR(raw_tau_vals)
    MALLOC_DLP_VAR(phase_deg_vals)
    MALLOC_DLP_VAR(raw_tau_threshold_vals)
    MALLOC_DLP_VAR(t_oet_spm_vals)
    MALLOC_DLP_VAR(t_ret_spm_vals)
    MALLOC_DLP_VAR(t_set_spm_vals)
    MALLOC_DLP_VAR(B_deg_vals)
    MALLOC_DLP_VAR(rho_dot_kms_vals)
    MALLOC_DLP_VAR(F_km_vals)
    MALLOC_DLP_VAR(D_km_vals)
    MALLOC_DLP_VAR(f_sky_hz_vals)
}
/*  End of rssringoccs_DLPCSV_Malloc.                                         */

/*  Undefine everything in case someone wants to #include this file.          */
#undef MALLOC_DLP_VAR
