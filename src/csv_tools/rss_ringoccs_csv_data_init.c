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
 *      Initialize all members of the CSV Data object to their zero values.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 4, 2024                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_CSVData typedef here, and function prototype given.           */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a Geo CSV to their default values.                   */
void rssringoccs_CSVData_Init(rssringoccs_CSVData *csv)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!csv)
        return;

    /*  All arrays should be empty, for now.                                  */
    csv->B_deg_vals = NULL;
    csv->D_km_vals = NULL;
    csv->f_sky_hz_vals = NULL;
    csv->p_norm_vals = NULL;
    csv->raw_tau_vals = NULL;
    csv->phase_deg_vals = NULL;
    csv->phi_deg_vals = NULL;
    csv->phi_rl_deg_vals = NULL;
    csv->raw_tau_threshold_vals = NULL;
    csv->rho_corr_pole_km_vals = NULL;
    csv->rho_corr_timing_km_vals = NULL;
    csv->rho_dot_kms_vals = NULL;
    csv->rho_km_vals = NULL;
    csv->rx_km_vals = NULL;
    csv->ry_km_vals = NULL;
    csv->rz_km_vals = NULL;
    csv->t_oet_spm_vals = NULL;
    csv->t_ret_spm_vals = NULL;
    csv->t_set_spm_vals = NULL;
    csv->tau_phase_deg_vals = NULL;
    csv->tau_power_vals = NULL;
    csv->tau_vals = NULL;
    csv->n_elements = 0;

    /*  The data entries should be nullified as well.                         */
    csv->geo = NULL;
    csv->cal = NULL;
    csv->dlp = NULL;
    csv->tau = NULL;

    /*  And no error has occurred yet.                                        */
    csv->error_message = NULL;
    csv->error_occurred = tmpl_False;

    /*  By default, do not use the deprecated file structure.                 */
    csv->use_deprecated = tmpl_False;

    /*  For chord occultations the geo file is incremented and decremented.   *
     *  By default these values should be zero.                               */
    csv->geo_increment = 0;
    csv->geo_decrement = 0;
}
/*  End of rssringoccs_CSVData_Init.                                          */
