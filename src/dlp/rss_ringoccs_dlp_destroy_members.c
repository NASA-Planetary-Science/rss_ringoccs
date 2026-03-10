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

/*  TMPL_FREE macro provided here, free's and nullifies a pointer.            */
#include <libtmpl/include/compat/tmpl_free.h>

/*  DLP object typedef found here.                                            */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  Function prototype / forward declaration.                                 */
extern void rssringoccs_DLP_Destroy_Members(rssringoccs_DLPObj * const dlp);

/*  Function for freeing all member of a dlp object except the error message. */
void rssringoccs_DLP_Destroy_Members(rssringoccs_DLPObj * const dlp)
{
    /*  If the input pointer is NULL, do not try to access it. Just return.   */
    if (!dlp)
        return;

    /*  Use the TMPL_FREE macro to free and Nullify all pointers.             */
    TMPL_FREE(dlp->rho_km_vals);
    TMPL_FREE(dlp->phi_deg_vals);
    TMPL_FREE(dlp->B_deg_vals);
    TMPL_FREE(dlp->D_km_vals);
    TMPL_FREE(dlp->f_sky_hz_vals);
    TMPL_FREE(dlp->rho_dot_kms_vals);
    TMPL_FREE(dlp->t_oet_spm_vals);
    TMPL_FREE(dlp->t_ret_spm_vals);
    TMPL_FREE(dlp->t_set_spm_vals);
    TMPL_FREE(dlp->rho_corr_pole_km_vals);
    TMPL_FREE(dlp->rho_corr_timing_km_vals);
    TMPL_FREE(dlp->phi_rl_deg_vals);
    TMPL_FREE(dlp->p_norm_vals);
    TMPL_FREE(dlp->phase_deg_vals);
    TMPL_FREE(dlp->raw_tau_threshold_vals);
    TMPL_FREE(dlp->rx_km_vals);
    TMPL_FREE(dlp->ry_km_vals);
    TMPL_FREE(dlp->rz_km_vals);
}
/*  End of rssringoccs_DLP_Destroy_Members.                                   */
