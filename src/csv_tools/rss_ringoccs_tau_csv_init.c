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
 *      Initialize all members of the Tau CSV object to their zero values.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_TauCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a Tau CSV to their default values.                   */
void rssringoccs_TauCSV_Init(rssringoccs_TauCSV *tau)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!tau)
        return;

    /*  All arrays should be empty, for now.                                  */
    tau->rho_km_vals = NULL;
    tau->rho_corr_pole_km_vals = NULL;
    tau->rho_corr_timing_km_vals = NULL;
    tau->phi_rl_deg_vals = NULL;
    tau->phi_ora_deg_vals = NULL;
    tau->power_vals = NULL;
    tau->tau_vals = NULL;
    tau->phase_deg_vals = NULL;
    tau->tau_threshold_vals = NULL;
    tau->t_oet_spm_vals = NULL;
    tau->t_ret_spm_vals = NULL;
    tau->t_set_spm_vals = NULL;
    tau->B_deg_vals = NULL;
    tau->error_message = NULL;
    tau->error_occurred = tmpl_False;
    tau->n_elements = 0;

    /*  And no error has occurred, yet.                                       */
    tau->error_message = NULL;
    tau->error_occurred = tmpl_False;

    /*  By default, do not use the deprecated file structure.                 */
    tau->use_deprecated = tmpl_False;
}
/*  End of rssringoccs_TauCSV_Init.                                           */
