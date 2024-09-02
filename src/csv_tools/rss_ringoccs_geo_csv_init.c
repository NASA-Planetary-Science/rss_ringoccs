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
 *      Initialize all members of the Geo CSV object to their zero values.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_GeoCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a Geo CSV to their default values.                   */
void rssringoccs_GeoCSV_Init(rssringoccs_GeoCSV *geo)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (geo == NULL)
        return;

    /*  All arrays should be empty, for now.                                  */
    geo->t_oet_spm_vals = NULL;
    geo->t_ret_spm_vals = NULL;
    geo->t_set_spm_vals = NULL;
    geo->rho_km_vals = NULL;
    geo->phi_rl_deg_vals = NULL;
    geo->phi_ora_deg_vals = NULL;
    geo->B_deg_vals = NULL;
    geo->D_km_vals = NULL;
    geo->rho_dot_kms_vals = NULL;
    geo->phi_rl_dot_kms_vals = NULL;
    geo->F_km_vals = NULL;
    geo->R_imp_km_vals = NULL;
    geo->rx_km_vals = NULL;
    geo->ry_km_vals = NULL;
    geo->rz_km_vals = NULL;
    geo->vx_kms_vals = NULL;
    geo->vy_kms_vals = NULL;
    geo->vz_kms_vals = NULL;
    geo->obs_spacecraft_lat_deg_vals = NULL;
    geo->n_elements = 0;

    /*  And no error has occurred yet.                                        */
    geo->error_message = NULL;
    geo->error_occurred = tmpl_False;

    /*  By default, do not use the deprecated file structure.                 */
    geo->use_deprecated = tmpl_False;
}
/*  End of rssringoccs_GeoCSV_Init.                                           */
