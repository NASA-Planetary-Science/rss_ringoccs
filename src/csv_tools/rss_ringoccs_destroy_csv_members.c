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
 *      Free all of the pointers in a CSVData object.                         *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 1, 2021                                               *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_CSVData typedef here, and function prototype given.           */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name has been defined.                                */
#ifdef DESTROY_CSV_VAR
#undef DESTROY_CSV_VAR
#endif

/*  Macro for freeing and nullifying the members of the CSV struct.           */
#define DESTROY_CSV_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of a rssringoccs_CSVData pointer except the            *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_Destroy_CSV_Members(rssringoccs_CSVData *csv)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (csv == NULL)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_CSV_VAR(csv->B_deg_vals)
    DESTROY_CSV_VAR(csv->D_km_vals)
    DESTROY_CSV_VAR(csv->f_sky_hz_vals)
    DESTROY_CSV_VAR(csv->p_norm_vals)
    DESTROY_CSV_VAR(csv->raw_tau_vals)
    DESTROY_CSV_VAR(csv->phase_deg_vals)
    DESTROY_CSV_VAR(csv->phi_deg_vals)
    DESTROY_CSV_VAR(csv->phi_rl_deg_vals)
    DESTROY_CSV_VAR(csv->raw_tau_threshold_vals)
    DESTROY_CSV_VAR(csv->rho_corr_pole_km_vals)
    DESTROY_CSV_VAR(csv->rho_corr_timing_km_vals)
    DESTROY_CSV_VAR(csv->rho_dot_kms_vals)
    DESTROY_CSV_VAR(csv->rho_km_vals)
    DESTROY_CSV_VAR(csv->rx_km_vals)
    DESTROY_CSV_VAR(csv->ry_km_vals)
    DESTROY_CSV_VAR(csv->rz_km_vals)
    DESTROY_CSV_VAR(csv->t_oet_spm_vals)
    DESTROY_CSV_VAR(csv->t_ret_spm_vals)
    DESTROY_CSV_VAR(csv->t_set_spm_vals)
    DESTROY_CSV_VAR(csv->tau_phase)
    DESTROY_CSV_VAR(csv->tau_power)
    DESTROY_CSV_VAR(csv->tau_vals)
}
/*  End of rssringoccs_Destroy_CSV_Members.                                   */

/*  Undefine the macro function.                                              */
#undef DESTROY_CSV_VAR
