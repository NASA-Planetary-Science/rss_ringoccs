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

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Header file with the DLP object definition.                               */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  NULL pointer is given here.                                               */
#include <stddef.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_DLP_Init(rssringoccs_DLPObj * const dlp);

/*  Function for initializing all members in a dlp object to NULL.            */
void rssringoccs_DLP_Init(rssringoccs_DLPObj * const dlp)
{
    /*  If the input is a NULL pointer there is nothing to be done. Abort.    */
    if (!dlp)
        return;

    /*  Initialize all pointers to NULL. This prevents things like double     *
     *  free's, leaking memory by calling malloc twice, etc. The functions    *
     *  that handle memory management with DLP objects will assume that       *
     *  either these pointers are NULL, or we successfully initialized using  *
     *  either malloc or calloc. Setting everything to NULL at the start      *
     *  helps reduce the chance of errors.                                    */
    dlp->rho_km_vals = NULL;
    dlp->phi_deg_vals = NULL;
    dlp->B_deg_vals = NULL;
    dlp->D_km_vals = NULL;
    dlp->f_sky_hz_vals = NULL;
    dlp->rho_dot_kms_vals = NULL;
    dlp->t_oet_spm_vals = NULL;
    dlp->t_ret_spm_vals = NULL;
    dlp->t_set_spm_vals = NULL;
    dlp->rho_corr_pole_km_vals = NULL;
    dlp->rho_corr_timing_km_vals = NULL;
    dlp->phi_rl_deg_vals = NULL;
    dlp->p_norm_vals = NULL;
    dlp->phase_deg_vals = NULL;
    dlp->raw_tau_threshold_vals = NULL;
    dlp->rx_km_vals = NULL;
    dlp->ry_km_vals = NULL;
    dlp->rz_km_vals = NULL;

    /*  The rho_km_vals pointer is currently NULL, set the displacement to 0. */
    dlp->dx_km = 0.0;

    /*  Set the indexing variables to be zero as well.                        */
    dlp->arr_size = 0;
    dlp->reference_count = 0;

    /*  Several functions will print messages throughout the computation if   *
     *  the verbose Boolean is set to True. Default is silent, set to False.  */
    dlp->verbose = tmpl_False;

    /*  Fresh DLP object, set the error variables to their zero values.       */
    dlp->error_occurred = tmpl_False;
    dlp->error_message = NULL;
}
/*  End of rssringoccs_DLP_Init.                                              */
