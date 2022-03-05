/******************************************************************************
 *                                 LICENSE                                    *
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
 *      Free all of the pointers in a TauCSV object.                          *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name is available.                                    */
#ifdef DESTROY_TAU_VAR
#undef DESTROY_TAU_VAR
#endif

/*  Macro for freeing and nullifying the members of the geo CSV structs.      */
#define DESTROY_TAU_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of a rssringoccs_TauCSV pointer except the             *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_Destroy_TauCSV_Members(rssringoccs_TauCSV *tau)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (tau == NULL)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_TAU_VAR(tau->rho_km_vals)
    DESTROY_TAU_VAR(tau->rho_corr_pole_km_vals)
    DESTROY_TAU_VAR(tau->rho_corr_timing_km_vals)
    DESTROY_TAU_VAR(tau->phi_rl_deg_vals)
    DESTROY_TAU_VAR(tau->phi_ora_deg_vals)
    DESTROY_TAU_VAR(tau->power_vals)
    DESTROY_TAU_VAR(tau->tau_vals)
    DESTROY_TAU_VAR(tau->phase_deg_vals)
    DESTROY_TAU_VAR(tau->tau_threshold_vals)
    DESTROY_TAU_VAR(tau->t_oet_spm_vals)
    DESTROY_TAU_VAR(tau->t_ret_spm_vals)
    DESTROY_TAU_VAR(tau->t_set_spm_vals)
    DESTROY_TAU_VAR(tau->B_deg_vals)
}
/*  End of rssringoccs_Destroy_TauCSV_Members.                                */

/*  Undefine the macro function.                                              */
#undef DESTROY_TAU_VAR
