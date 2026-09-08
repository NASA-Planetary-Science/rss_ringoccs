/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *      Free all of the pointers in a MergedCSVData object.                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 30, 2024                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  TMPL_FREE macro found here.                                               */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_MergedCSVData typedef here, and function prototype given.     */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Free's all members of a rssringoccs_MergedCSVData pointer except the      *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_MergedCSVData_Destroy_Members(rssringoccs_MergedCSVData *dlpm)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlpm)
        return;

    /*  Destroy every variable except the error_message.                      */
    TMPL_FREE(dlpm->B_deg_vals);
    TMPL_FREE(dlpm->D_km_vals);
    TMPL_FREE(dlpm->f_sky_hz_vals);
    TMPL_FREE(dlpm->p_norm_vals);
    TMPL_FREE(dlpm->raw_tau_vals);
    TMPL_FREE(dlpm->phase_deg_vals);
    TMPL_FREE(dlpm->phi_deg_vals);
    TMPL_FREE(dlpm->phi_rl_deg_vals);
    TMPL_FREE(dlpm->raw_tau_threshold_vals);
    TMPL_FREE(dlpm->rho_corr_pole_km_vals);
    TMPL_FREE(dlpm->rho_corr_timing_km_vals);
    TMPL_FREE(dlpm->rho_dot_kms_vals);
    TMPL_FREE(dlpm->rho_km_vals);
    TMPL_FREE(dlpm->rx_km_vals);
    TMPL_FREE(dlpm->ry_km_vals);
    TMPL_FREE(dlpm->rz_km_vals);
    TMPL_FREE(dlpm->t_oet_spm_vals);
    TMPL_FREE(dlpm->t_ret_spm_vals);
    TMPL_FREE(dlpm->t_set_spm_vals);
    TMPL_FREE(dlpm->tau_phase_deg_vals);
    TMPL_FREE(dlpm->tau_power_vals);
    TMPL_FREE(dlpm->tau_vals);
}
/*  End of rssringoccs_MergedCSVData_Destroy_Members.                         */
