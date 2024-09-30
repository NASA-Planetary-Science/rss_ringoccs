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
 *      Initialize all members of the DLPM CSV object to their zero values.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 30, 2024                                            *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_MergedCSVData typedef here, and function prototype given.     */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a DLPM CSV object to their default values.           */
void rssringoccs_MergedCSVData_Init(rssringoccs_MergedCSVData *dlpm)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlpm)
        return;

    /*  All arrays should be empty, for now.                                  */
    dlpm->B_deg_vals = NULL;
    dlpm->D_km_vals = NULL;
    dlpm->f_sky_hz_vals = NULL;
    dlpm->p_norm_vals = NULL;
    dlpm->raw_tau_vals = NULL;
    dlpm->phase_deg_vals = NULL;
    dlpm->phi_deg_vals = NULL;
    dlpm->phi_rl_deg_vals = NULL;
    dlpm->raw_tau_threshold_vals = NULL;
    dlpm->rho_corr_pole_km_vals = NULL;
    dlpm->rho_corr_timing_km_vals = NULL;
    dlpm->rho_dot_kms_vals = NULL;
    dlpm->rho_km_vals = NULL;
    dlpm->rx_km_vals = NULL;
    dlpm->ry_km_vals = NULL;
    dlpm->rz_km_vals = NULL;
    dlpm->t_oet_spm_vals = NULL;
    dlpm->t_ret_spm_vals = NULL;
    dlpm->t_set_spm_vals = NULL;
    dlpm->tau_phase_deg_vals = NULL;
    dlpm->tau_power_vals = NULL;
    dlpm->tau_vals = NULL;
    dlpm->n_elements = 0;

    /*  And no error has occurred, yet.                                       */
    dlpm->error_message = NULL;
    dlpm->error_occurred = tmpl_False;

    /*  No history for the object, either.                                    */
    dlpm->history = NULL;
}
/*  End of rssringoccs_MergedCSVData_Init.                                    */
