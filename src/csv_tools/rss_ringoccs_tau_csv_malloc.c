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
 *      Function for allocating memory for a Tau CSV object.                  *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  malloc is found here.                                                     */
#include <stdlib.h>

/*  libtmpl provided Booleans and line count.                                 */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_TauCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if the macro name is available.                                     */
#ifdef MALLOC_TAU_VAR
#undef MALLOC_TAU_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_TAU_VAR(var)                                                    \
    tau->var = malloc(sizeof(*tau->var) * tau->n_elements);                    \
    if (!tau->var)                                                             \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message =                                                   \
            "\nError Encountered: rss_ringoccs\n"                              \
            "\trssringoccs_TauCSV_Malloc\n\n"                                  \
            "malloc returned NULL. Failed to allocate memory for " #var ".\n"; \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_TauCSV_Destroy_Members(tau);                               \
        return;                                                                \
    }

/*  Function for allocating memory to a Tau CSV based on a CSV file pointer.  */
void rssringoccs_TauCSV_Malloc(rssringoccs_TauCSV *tau, FILE *fp)
{
    /*  If the input Tau object is NULL there is nothing to do. Return.       */
    if (!tau)
        return;

    /*  Similarly if an error already occurred. Abort the computation.        */
    if (tau->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_TauCSV_Malloc\n\n"
            "Input file is NULL.\n";

        return;
    }

    /*  Count the number of lines in the CSV.                                 */
    tau->n_elements = tmpl_Line_Count(fp);

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (tau->n_elements == (size_t)0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_TauCSV_Malloc\n\n"
            "n_elements is zero, nothing to malloc.\n";

        return;
    }

    /*  Use the MALLOC_TAU_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_TAU_VAR(rho_km_vals)
    MALLOC_TAU_VAR(rho_corr_pole_km_vals)
    MALLOC_TAU_VAR(rho_corr_timing_km_vals)
    MALLOC_TAU_VAR(phi_rl_deg_vals)
    MALLOC_TAU_VAR(phi_ora_deg_vals)
    MALLOC_TAU_VAR(tau_vals)
    MALLOC_TAU_VAR(phase_deg_vals)
    MALLOC_TAU_VAR(tau_threshold_vals)
    MALLOC_TAU_VAR(t_oet_spm_vals)
    MALLOC_TAU_VAR(t_ret_spm_vals)
    MALLOC_TAU_VAR(t_set_spm_vals)
    MALLOC_TAU_VAR(B_deg_vals)

    /*  If we're using the older deprecated format, there are 13 columns.     *
     *  Allocate memory for power_vals as well.                               */
    if (!tau->use_deprecated)
    {
        MALLOC_TAU_VAR(power_vals)
    }
}
/*  End of rssringoccs_TauCSV_Malloc.                                         */

/*  Undefine everything in case someone wants to #include this file.          */
#undef MALLOC_TAU_VAR
