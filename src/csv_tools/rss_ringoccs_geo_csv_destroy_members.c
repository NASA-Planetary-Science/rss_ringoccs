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
 *      Free all of the pointers in a GeoCSV object.                          *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Macro for freeing a pointer and setting it to NULL.                       */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_GeoCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Free's all members of a rssringoccs_GeoCSV pointer except the             *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_GeoCSV_Destroy_Members(rssringoccs_GeoCSV * const geo)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!geo)
        return;

    /*  Destroy every variable except the error_message.                      */
    TMPL_FREE(geo->t_oet_spm_vals);
    TMPL_FREE(geo->t_ret_spm_vals);
    TMPL_FREE(geo->t_set_spm_vals);
    TMPL_FREE(geo->rho_km_vals);
    TMPL_FREE(geo->phi_rl_deg_vals);
    TMPL_FREE(geo->phi_ora_deg_vals);
    TMPL_FREE(geo->B_deg_vals);
    TMPL_FREE(geo->D_km_vals);
    TMPL_FREE(geo->rho_dot_kms_vals);
    TMPL_FREE(geo->phi_rl_dot_kms_vals);
    TMPL_FREE(geo->F_km_vals);
    TMPL_FREE(geo->R_imp_km_vals);
    TMPL_FREE(geo->rx_km_vals);
    TMPL_FREE(geo->ry_km_vals);
    TMPL_FREE(geo->rz_km_vals);
    TMPL_FREE(geo->vx_kms_vals);
    TMPL_FREE(geo->vy_kms_vals);
    TMPL_FREE(geo->vz_kms_vals);
    TMPL_FREE(geo->obs_spacecraft_lat_deg_vals);
    TMPL_FREE(geo->history);
}
/*  End of rssringoccs_GeoCSV_Destroy_Members.                                */
