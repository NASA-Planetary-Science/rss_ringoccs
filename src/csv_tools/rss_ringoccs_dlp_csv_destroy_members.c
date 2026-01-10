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
 *      Free all of the pointers in a DLPCSV object.                          *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Macro for freeing a pointer and setting it to NULL.                       */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Free's all members of a rssringoccs_DLPCSV pointer except the             *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_DLPCSV_Destroy_Members(rssringoccs_DLPCSV * const dlp)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlp)
        return;

    /*  Destroy every variable except the error_message.                      */
    TMPL_FREE(dlp->rho_km_vals);
    TMPL_FREE(dlp->rho_corr_pole_km_vals);
    TMPL_FREE(dlp->rho_corr_timing_km_vals);
    TMPL_FREE(dlp->phi_rl_deg_vals);
    TMPL_FREE(dlp->phi_ora_deg_vals);
    TMPL_FREE(dlp->p_norm_vals);
    TMPL_FREE(dlp->raw_tau_vals);
    TMPL_FREE(dlp->phase_deg_vals);
    TMPL_FREE(dlp->raw_tau_threshold_vals);
    TMPL_FREE(dlp->t_oet_spm_vals);
    TMPL_FREE(dlp->t_ret_spm_vals);
    TMPL_FREE(dlp->t_set_spm_vals);
    TMPL_FREE(dlp->B_deg_vals);
}
/*  End of rssringoccs_DLPCSV_Destroy_Members.                                */
