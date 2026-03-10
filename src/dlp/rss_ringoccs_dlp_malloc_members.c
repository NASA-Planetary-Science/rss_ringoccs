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
 *  Author:     Ryan Maguire                                                  *
 *  Date:       February 27, 2026                                             *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  TMPL_MALLOC macro found here, providing C vs. C++ compatibility.          */
#include <libtmpl/include/compat/tmpl_malloc.h>

/*  Header file with the DLP object definition.                               */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_DLP_Malloc_Members(rssringoccs_DLPObj * const dlp);

/*  Use this macro to save on repetitive code. It checks if dlp->var is NULL, *
 *  attempts to malloc memory for dlp->var if it is, and then checks to see   *
 *  if malloc failed.                                                         */
#define MALLOC_DLP_VAR(var)                                                    \
    do {                                                                       \
        /*  Check if the variable is not NULL. It should be at the start.    */\
        if (dlp->var)                                                          \
        {                                                                      \
            dlp->error_occurred = tmpl_True;                                   \
            dlp->error_message =                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\trssringoccs_DLP_Malloc_Members\n\n"                       \
                "\r"#var" is not NULL. It is likely you've already set the\n"  \
                "\rdata for this DLP object.\n\n";                             \
            return;                                                            \
        }                                                                      \
                                                                               \
        /*  Allocate memory for the variable.                                */\
        dlp->var = TMPL_MALLOC(double, dlp->arr_size);                         \
                                                                               \
        /*  Check if malloc failed.                                          */\
        if (!dlp->var)                                                         \
        {                                                                      \
            dlp->error_occurred = tmpl_True;                                   \
            dlp->error_message =                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\trssringoccs_DLP_Malloc_Members\n\n"                       \
                "\rMalloc failed and returned NULL for "#var".\n\n";           \
            return;                                                            \
        }                                                                      \
    } while (0)
/*  End of the MALLOC_DLP_VAR macro.                                          */

/*  Function for allocating memory for all of the dlp variables.              */
void rssringoccs_DLP_Malloc_Members(rssringoccs_DLPObj * const dlp)
{
    if (!dlp)
        return;

    if (dlp->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (dlp->verbose)
        puts("\r\tDLP: Allocating memory for core arrays...");

    if (dlp->arr_size == 0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Malloc_Members\n\n"
            "\rInput dlp has arr_size = 0, nothing to allocate.\n\n";

        return;
    }

    /*  Allocate memory for each of the arrays in a DLP object.               */
    MALLOC_DLP_VAR(rho_km_vals);
    MALLOC_DLP_VAR(phi_deg_vals);
    MALLOC_DLP_VAR(B_deg_vals);
    MALLOC_DLP_VAR(D_km_vals);
    MALLOC_DLP_VAR(f_sky_hz_vals);
    MALLOC_DLP_VAR(rho_dot_kms_vals);
    MALLOC_DLP_VAR(t_oet_spm_vals);
    MALLOC_DLP_VAR(t_ret_spm_vals);
    MALLOC_DLP_VAR(t_set_spm_vals);
    MALLOC_DLP_VAR(rho_corr_pole_km_vals);
    MALLOC_DLP_VAR(rho_corr_timing_km_vals);
    MALLOC_DLP_VAR(phi_rl_deg_vals);
    MALLOC_DLP_VAR(p_norm_vals);
    MALLOC_DLP_VAR(phase_deg_vals);
    MALLOC_DLP_VAR(raw_tau_threshold_vals);
    MALLOC_DLP_VAR(rx_km_vals);
    MALLOC_DLP_VAR(ry_km_vals);
    MALLOC_DLP_VAR(rz_km_vals);
}
/*  End of rssringoccs_DLP_Malloc_Members.                                    */

/*  Undefine this in case someone wants to #include this file.                */
#undef MALLOC_DLP_VAR
