/******************************************************************************
 *                                 LICENSE                                    *
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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 5, 2020                                               *
 ******************************************************************************/

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Macro for freeing and nullifying the members of the geo CSV structs.      */
#define DESTROY_TAU_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Function for freeing all member of a tau object except the error message. */
void rssringoccs_Destroy_Tau_Members(rssringoccs_TAUObj *tau)
{
    /*  If the input pointer is NULL, do not try to access it. Just return.   */
    if (tau == NULL)
        return;

    /*  Use the DESTROY_TAU_VAR macro to free and NULLify all pointers.       */
    DESTROY_TAU_VAR(tau->rho_km_vals);
    DESTROY_TAU_VAR(tau->F_km_vals);
    DESTROY_TAU_VAR(tau->phi_rad_vals);
    DESTROY_TAU_VAR(tau->k_vals);
    DESTROY_TAU_VAR(tau->f_sky_hz_vals);
    DESTROY_TAU_VAR(tau->rho_dot_kms_vals);
    DESTROY_TAU_VAR(tau->raw_tau_threshold_vals);
    DESTROY_TAU_VAR(tau->B_rad_vals);
    DESTROY_TAU_VAR(tau->D_km_vals);
    DESTROY_TAU_VAR(tau->w_km_vals);
    DESTROY_TAU_VAR(tau->t_oet_spm_vals);
    DESTROY_TAU_VAR(tau->t_ret_spm_vals);
    DESTROY_TAU_VAR(tau->t_set_spm_vals);
    DESTROY_TAU_VAR(tau->rho_corr_pole_km_vals);
    DESTROY_TAU_VAR(tau->rho_corr_timing_km_vals);
    DESTROY_TAU_VAR(tau->tau_threshold_vals);
    DESTROY_TAU_VAR(tau->phi_rl_rad_vals);
    DESTROY_TAU_VAR(tau->p_norm_vals);
    DESTROY_TAU_VAR(tau->p_norm_fwd_vals);
    DESTROY_TAU_VAR(tau->power_vals);
    DESTROY_TAU_VAR(tau->phase_rad_vals);
    DESTROY_TAU_VAR(tau->phase_fwd_vals);
    DESTROY_TAU_VAR(tau->phase_vals);
    DESTROY_TAU_VAR(tau->tau_vals);
    DESTROY_TAU_VAR(tau->tau_fwd_vals);
    DESTROY_TAU_VAR(tau->rx_km_vals);
    DESTROY_TAU_VAR(tau->ry_km_vals);
    DESTROY_TAU_VAR(tau->rz_km_vals);
    DESTROY_TAU_VAR(tau->wtype);
    DESTROY_TAU_VAR(tau->psitype);
    DESTROY_TAU_VAR(tau->T_in);
    DESTROY_TAU_VAR(tau->T_out);
    DESTROY_TAU_VAR(tau->T_fwd);
}
/*  End of rssringoccs_Destroy_Tau_Members.                                   */
