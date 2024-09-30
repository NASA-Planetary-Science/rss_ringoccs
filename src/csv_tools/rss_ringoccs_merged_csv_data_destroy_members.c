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

/*  rssringoccs_MergedCSVData typedef here, and function prototype given.     */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name is available.                                    */
#ifdef DESTROY_DLPM_VAR
#undef DESTROY_DLPM_VAR
#endif

/*  Macro for freeing and nullifying the members of the DLPM CSV struct.      */
#define DESTROY_DLPM_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of a rssringoccs_MergedCSVData pointer except the      *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_MergedCSVData_Destroy_Members(rssringoccs_MergedCSVData *dlpm)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlpm)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_DLPM_VAR(dlpm->B_deg_vals);
    DESTROY_DLPM_VAR(dlpm->D_km_vals);
    DESTROY_DLPM_VAR(dlpm->f_sky_hz_vals);
    DESTROY_DLPM_VAR(dlpm->p_norm_vals);
    DESTROY_DLPM_VAR(dlpm->raw_tau_vals);
    DESTROY_DLPM_VAR(dlpm->phase_deg_vals);
    DESTROY_DLPM_VAR(dlpm->phi_deg_vals);
    DESTROY_DLPM_VAR(dlpm->phi_rl_deg_vals);
    DESTROY_DLPM_VAR(dlpm->raw_tau_threshold_vals);
    DESTROY_DLPM_VAR(dlpm->rho_corr_pole_km_vals);
    DESTROY_DLPM_VAR(dlpm->rho_corr_timing_km_vals);
    DESTROY_DLPM_VAR(dlpm->rho_dot_kms_vals);
    DESTROY_DLPM_VAR(dlpm->rho_km_vals);
    DESTROY_DLPM_VAR(dlpm->rx_km_vals);
    DESTROY_DLPM_VAR(dlpm->ry_km_vals);
    DESTROY_DLPM_VAR(dlpm->rz_km_vals);
    DESTROY_DLPM_VAR(dlpm->t_oet_spm_vals);
    DESTROY_DLPM_VAR(dlpm->t_ret_spm_vals);
    DESTROY_DLPM_VAR(dlpm->t_set_spm_vals);
    DESTROY_DLPM_VAR(dlpm->tau_phase_deg_vals);
    DESTROY_DLPM_VAR(dlpm->tau_power_vals);
    DESTROY_DLPM_VAR(dlpm->tau_vals);
}
/*  End of rssringoccs_MergedCSVData_Destroy_Members.                         */

/*  Undefine everything in case someone wants to #include this file.          */
#undef DESTROY_DLPM_VAR
