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
 *      Function for allocating memory for a DLPM CSV object.                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 30, 2024                                            *
 ******************************************************************************/

/*  malloc is found here.                                                     */
#include <stdlib.h>

/*  libtmpl provided Booleans, string duplicate, and line count.              */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_MergedCSVData typedef here, and function prototype given.     */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if the macro name is available.                                     */
#ifdef MALLOC_DLPM_VAR
#undef MALLOC_DLPM_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_DLPM_VAR(var)                                                   \
    var = malloc(sizeof(*var) * dlpm->n_elements);                             \
    if (!var)                                                                  \
    {                                                                          \
        dlpm->error_occurred = tmpl_True;                                      \
        dlpm->error_message = tmpl_String_Duplicate(                           \
            "Error Encountered: rss_ringoccs\n"                                \
            "\trssringoccs_MergedCSVData_Malloc\n\n"                           \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_MergedCSVData_Destroy_Members(dlpm);                       \
        return;                                                                \
    }

/*  Function for allocating memory to a DLPM CSV based on a CSV file pointer. */
void rssringoccs_MergedCSVData_Malloc(rssringoccs_MergedCSVData *dlpm, FILE *fp)
{
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
            "\trssringoccs_MergedCSVData_Malloc\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  Count the number of lines in the CSV.                                 */
    dlpm->n_elements = tmpl_Line_Count(fp);

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (dlpm->n_elements == (size_t)0)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Malloc\n\n"
            "n_elements is zero, nothing to malloc. Aborting.\n"
        );

        return;
    }

    /*  Use the MALLOC_DLPM_VAR macro function to allocate memory and check   *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_DLPM_VAR(dlpm->B_deg_vals)
    MALLOC_DLPM_VAR(dlpm->D_km_vals)
    MALLOC_DLPM_VAR(dlpm->f_sky_hz_vals)
    MALLOC_DLPM_VAR(dlpm->p_norm_vals)
    MALLOC_DLPM_VAR(dlpm->raw_tau_vals)
    MALLOC_DLPM_VAR(dlpm->phase_deg_vals)
    MALLOC_DLPM_VAR(dlpm->phi_deg_vals)
    MALLOC_DLPM_VAR(dlpm->phi_rl_deg_vals)
    MALLOC_DLPM_VAR(dlpm->raw_tau_threshold_vals)
    MALLOC_DLPM_VAR(dlpm->rho_corr_pole_km_vals)
    MALLOC_DLPM_VAR(dlpm->rho_corr_timing_km_vals)
    MALLOC_DLPM_VAR(dlpm->rho_dot_kms_vals)
    MALLOC_DLPM_VAR(dlpm->rho_km_vals)
    MALLOC_DLPM_VAR(dlpm->rx_km_vals)
    MALLOC_DLPM_VAR(dlpm->ry_km_vals)
    MALLOC_DLPM_VAR(dlpm->rz_km_vals)
    MALLOC_DLPM_VAR(dlpm->t_oet_spm_vals)
    MALLOC_DLPM_VAR(dlpm->t_ret_spm_vals)
    MALLOC_DLPM_VAR(dlpm->t_set_spm_vals)
    MALLOC_DLPM_VAR(dlpm->tau_phase_deg_vals)
    MALLOC_DLPM_VAR(dlpm->tau_power_vals)
    MALLOC_DLPM_VAR(dlpm->tau_vals)
}
/*  End of rssringoccs_MergedCSVData_Malloc.                                  */

/*  Undefine everything in case someone wants to #include this file.          */
#undef MALLOC_DLPM_VAR
