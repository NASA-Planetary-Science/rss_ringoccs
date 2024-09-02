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
 *      Function for allocating memory for a Cal CSV object.                  *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  malloc is found here.                                                     */
#include <stdlib.h>

/*  libtmpl provides Booleans, string duplicate, and line count.              */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_CalCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name has been defined. Clear it if so, as we'll be    *
 *  using this macro to malloc the variables in the CAL struct.               */
#ifdef MALLOC_CAL_VAR
#undef MALLOC_CAL_VAR
#endif

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_CAL_VAR(var)                                                    \
    cal->var = malloc(sizeof(*cal->var) * cal->n_elements);                    \
    if (!cal->var)                                                             \
    {                                                                          \
        cal->error_occurred = tmpl_True;                                       \
        cal->error_message = tmpl_strdup(                                      \
            "Error Encountered: rss_ringoccs\n"                                \
            "\trssringoccs_CalCSV_Malloc\n\n"                                  \
            "Malloc returned NULL. Failed to allocate memory for " #var ".\n"  \
            "Aborting computation and returning.\n"                            \
        );                                                                     \
                                                                               \
        /*  Free the variables that have been malloc'd so far.               */\
        rssringoccs_CalCSV_Destroy_Members(cal);                               \
        return;                                                                \
    }

/*  Function of allocating memory to a Cal CSV based on a CSV file pointer.   */
void rssringoccs_CalCSV_Malloc(rssringoccs_CalCSV *cal, FILE *fp)
{
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
            "\trssringoccs_CalCSV_Malloc\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  Count the number of lines in the CSV.                                 */
    cal->n_elements = tmpl_Line_Count(fp);

    /*  There needs to be at least one row in the CSV file. If not, treat     *
     *  this as an error. It is likely the file is corrupted.                 */
    if (cal->n_elements == (size_t)0)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Malloc\n\n"
            "n_elements is zero, nothing to malloc. Aborting.\n"
        );

        return;
    }

    /*  Use the MALLOC_CAL_VAR macro function to allocate memory and check    *
     *  for errors. This macro ends with an if-then statement, and ends in    *
     *  curly braces {}, hence no need for a semi-colon here.                 */
    MALLOC_CAL_VAR(t_oet_spm_vals)
    MALLOC_CAL_VAR(f_sky_pred_vals)
    MALLOC_CAL_VAR(f_sky_resid_fit_vals)
    MALLOC_CAL_VAR(p_free_vals)
}
/*  End of rssringoccs_CalCSV_Malloc.                                         */

/*  Undefine everything in case someone wants to #include this file.          */
#undef MALLOC_CAL_VAR
