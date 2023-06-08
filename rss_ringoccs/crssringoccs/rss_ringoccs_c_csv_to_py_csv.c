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
 ******************************************************************************/
#include "crss_ringoccs.h"

void rssringoccs_C_CSV_to_Py_CSV(PyCSVObj *py_csv, rssringoccs_CSVData *csv)
{
    if (csv == NULL)
        return;

    if (csv->error_occurred)
        return;

    if (py_csv == NULL)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_CSV_to_Py_CSV\n\n"
            "\rInput py_csv is NULL. Aborting.n"
        );
        return;
    }

    set_var(&py_csv->rho_km_vals, &csv->rho_km_vals, csv->n_elements);
    set_var(&py_csv->B_rad_vals, &csv->B_rad_vals, csv->n_elements);
    set_var(&py_csv->D_km_vals, &csv->D_km_vals, csv->n_elements);
    set_var(&py_csv->f_sky_hz_vals, &csv->f_sky_hz_vals, csv->n_elements);
    set_var(&py_csv->p_norm_vals, &csv->p_norm_vals, csv->n_elements);
    set_var(&py_csv->phase_rad_vals, &csv->phase_rad_vals, csv->n_elements);
    set_var(&py_csv->phi_rad_vals, &csv->phi_rad_vals, csv->n_elements);
    set_var(&py_csv->phi_rl_rad_vals, &csv->phi_rl_rad_vals, csv->n_elements);
    set_var(&py_csv->rho_dot_kms_vals, &csv->rho_dot_kms_vals, csv->n_elements);
    set_var(&py_csv->t_oet_spm_vals, &csv->t_oet_spm_vals, csv->n_elements);
    set_var(&py_csv->t_ret_spm_vals, &csv->t_ret_spm_vals, csv->n_elements);
    set_var(&py_csv->t_set_spm_vals, &csv->t_set_spm_vals, csv->n_elements);
    set_var(&py_csv->tau_vals, &csv->tau_vals, csv->n_elements);
    set_var(&py_csv->rx_km_vals, &csv->rx_km_vals, csv->n_elements);
    set_var(&py_csv->ry_km_vals, &csv->ry_km_vals, csv->n_elements);
    set_var(&py_csv->rz_km_vals, &csv->rz_km_vals, csv->n_elements);
    set_var(&py_csv->raw_tau_threshold_vals, &csv->raw_tau_threshold_vals, csv->n_elements);
    set_var(&py_csv->rho_corr_pole_km_vals, &csv->rho_corr_pole_km_vals, csv->n_elements);
    set_var(&py_csv->rho_corr_timing_km_vals, &csv->rho_corr_timing_km_vals, csv->n_elements);
}
