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
 *      Free all of the pointers in a Uranus DLPCSV object.                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 24, 2024                                            *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_UranusDLPCSV typedef here, and function prototype given.      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name is available.                                    */
#ifdef DESTROY_DLP_VAR
#undef DESTROY_DLP_VAR
#endif

/*  Macro for freeing and nullifying the members of the dlp CSV struct.       */
#define DESTROY_DLP_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of a rssringoccs_UranusDLPCSV pointer except the       *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_UranusDLPCSV_Destroy_Members(rssringoccs_UranusDLPCSV *dlp)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlp)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_DLP_VAR(dlp->rho_km_vals)
    DESTROY_DLP_VAR(dlp->rho_corr_pole_km_vals)
    DESTROY_DLP_VAR(dlp->rho_corr_timing_km_vals)
    DESTROY_DLP_VAR(dlp->phi_rl_deg_vals)
    DESTROY_DLP_VAR(dlp->phi_deg_vals)
    DESTROY_DLP_VAR(dlp->p_norm_vals)
    DESTROY_DLP_VAR(dlp->raw_tau_vals)
    DESTROY_DLP_VAR(dlp->phase_deg_vals)
    DESTROY_DLP_VAR(dlp->raw_tau_threshold_vals)
    DESTROY_DLP_VAR(dlp->t_oet_spm_vals)
    DESTROY_DLP_VAR(dlp->t_ret_spm_vals)
    DESTROY_DLP_VAR(dlp->t_set_spm_vals)
    DESTROY_DLP_VAR(dlp->B_deg_vals)
    DESTROY_DLP_VAR(dlp->rho_dot_kms_vals)
    DESTROY_DLP_VAR(dlp->F_km_vals)
    DESTROY_DLP_VAR(dlp->D_km_vals)
    DESTROY_DLP_VAR(dlp->f_sky_hz_vals)
}
/*  End of rssringoccs_UranusDLPCSV_Destroy_Members.                          */

/*  Undefine everything in case someone wants to #include this file.          */
#undef DESTROY_DLP_VAR
