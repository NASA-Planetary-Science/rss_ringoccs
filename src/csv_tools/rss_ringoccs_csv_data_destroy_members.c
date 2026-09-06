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

/*  TMPL_FREE macro found here.                                               */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_CSVData typedef here, and function prototype given.           */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Free's all members of a rssringoccs_CSVData pointer except the            *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_CSVData_Destroy_Members(rssringoccs_CSVData *csv)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!csv)
        return;

    /*  Destroy every variable except the error_message.                      */
    TMPL_FREE(csv->B_deg_vals);
    TMPL_FREE(csv->D_km_vals);
    TMPL_FREE(csv->f_sky_hz_vals);
    TMPL_FREE(csv->p_norm_vals);
    TMPL_FREE(csv->raw_tau_vals);
    TMPL_FREE(csv->phase_deg_vals);
    TMPL_FREE(csv->phi_deg_vals);
    TMPL_FREE(csv->phi_rl_deg_vals);
    TMPL_FREE(csv->raw_tau_threshold_vals);
    TMPL_FREE(csv->rho_corr_pole_km_vals);
    TMPL_FREE(csv->rho_corr_timing_km_vals);
    TMPL_FREE(csv->rho_dot_kms_vals);
    TMPL_FREE(csv->rho_km_vals);
    TMPL_FREE(csv->rx_km_vals);
    TMPL_FREE(csv->ry_km_vals);
    TMPL_FREE(csv->rz_km_vals);
    TMPL_FREE(csv->t_oet_spm_vals);
    TMPL_FREE(csv->t_ret_spm_vals);
    TMPL_FREE(csv->t_set_spm_vals);
    TMPL_FREE(csv->tau_phase_deg_vals);
    TMPL_FREE(csv->tau_power_vals);
    TMPL_FREE(csv->tau_vals);

    /*  Destroy the CSV data if they exist.                                   */
    rssringoccs_GeoCSV_Destroy(&(csv->geo));
    rssringoccs_CalCSV_Destroy(&(csv->cal));
    rssringoccs_DLPCSV_Destroy(&(csv->dlp));
    rssringoccs_TauCSV_Destroy(&(csv->tau));
}
/*  End of rssringoccs_Destroy_CSV_Members.                                   */
