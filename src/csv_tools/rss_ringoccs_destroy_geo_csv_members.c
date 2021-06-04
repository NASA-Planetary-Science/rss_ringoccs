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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name is available.                                    */
#ifdef DESTROY_GEO_VAR
#undef DESTROY_GEO_VAR
#endif

/*  Macro for freeing and nullifying the members of the geo CSV structs.      */
#define DESTROY_GEO_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of a rssringoccs_GeoCSV pointer except the             *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_Destroy_GeoCSV_Members(rssringoccs_GeoCSV *geo)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (geo == NULL)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_GEO_VAR(geo->t_oet_spm_vals);
    DESTROY_GEO_VAR(geo->t_ret_spm_vals);
    DESTROY_GEO_VAR(geo->t_set_spm_vals);
    DESTROY_GEO_VAR(geo->rho_km_vals);
    DESTROY_GEO_VAR(geo->phi_rl_deg_vals);
    DESTROY_GEO_VAR(geo->phi_ora_deg_vals);
    DESTROY_GEO_VAR(geo->B_deg_vals);
    DESTROY_GEO_VAR(geo->D_km_vals);
    DESTROY_GEO_VAR(geo->rho_dot_kms_vals);
    DESTROY_GEO_VAR(geo->phi_rl_dot_kms_vals);
    DESTROY_GEO_VAR(geo->F_km_vals);
    DESTROY_GEO_VAR(geo->R_imp_km_vals);
    DESTROY_GEO_VAR(geo->rx_km_vals);
    DESTROY_GEO_VAR(geo->ry_km_vals);
    DESTROY_GEO_VAR(geo->rz_km_vals);
    DESTROY_GEO_VAR(geo->vx_kms_vals);
    DESTROY_GEO_VAR(geo->vy_kms_vals);
    DESTROY_GEO_VAR(geo->vz_kms_vals);
    DESTROY_GEO_VAR(geo->obs_spacecract_lat_deg_vals);
}
/*  End of rssringoccs_Destroy_GeoCSV_Members.                                */

/*  Undefine the macro function.                                              */
#undef DESTROY_GEO_VAR
