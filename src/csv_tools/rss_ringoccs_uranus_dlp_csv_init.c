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
 *      Initialize all members of the DLP CSV object to their zero values.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 24, 2024                                            *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_UranusDLPCSV typedef here, and function prototype given.      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a DLP CSV to their default values.                   */
void rssringoccs_UranusDLPCSV_Init(rssringoccs_UranusDLPCSV *dlp)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!dlp)
        return;

    /*  All arrays should be empty, for now.                                  */
    dlp->rho_km_vals = NULL;
    dlp->rho_corr_pole_km_vals = NULL;
    dlp->rho_corr_timing_km_vals = NULL;
    dlp->phi_rl_deg_vals = NULL;
    dlp->phi_deg_vals = NULL;
    dlp->p_norm_vals = NULL;
    dlp->raw_tau_vals = NULL;
    dlp->phase_deg_vals = NULL;
    dlp->raw_tau_threshold_vals = NULL;
    dlp->t_oet_spm_vals = NULL;
    dlp->t_ret_spm_vals = NULL;
    dlp->t_set_spm_vals = NULL;
    dlp->B_deg_vals = NULL;
    dlp->rho_dot_kms_vals = NULL;
    dlp->F_km_vals = NULL;
    dlp->D_km_vals = NULL;
    dlp->f_sky_hz_vals = NULL;
    dlp->n_elements = 0;

    /*  And no error has occurred, yet.                                       */
    dlp->error_message = NULL;
    dlp->error_occurred = tmpl_False;

    /*  Most DLP data sets for Uranus are in radians.                         */
    dlp->in_radians = tmpl_True;
}
/*  End of rssringoccs_UranusDLPCSV_Init.                                     */
