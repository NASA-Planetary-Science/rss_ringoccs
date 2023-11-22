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
 ******************************************************************************/

/*  NULL is defined here.                                                     */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_strdup function declared here.                                       */
#include <libtmpl/include/tmpl_string.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Macro for the crssringoccs_set_var function to shorten the syntax.        */
#define SET_VAR(a) crssringoccs_Set_Var(&py_csv->a, csv->a, csv->n_elements)

/*  Converts a C CSV struct to a Python CSV object.                           */
void crssringoccs_C_CSV_to_Py_CSV(PyCSVObj *py_csv, rssringoccs_CSVData *csv)
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
            "\r\tcrssringoccs_C_CSV_to_Py_CSV\n\n"
            "\rInput py_csv is NULL. Aborting.n"
        );
        return;
    }

    SET_VAR(rho_km_vals);
    SET_VAR(B_deg_vals);
    SET_VAR(D_km_vals);
    SET_VAR(f_sky_hz_vals);
    SET_VAR(p_norm_vals);
    SET_VAR(raw_tau_vals);
    SET_VAR(phase_deg_vals);
    SET_VAR(phi_deg_vals);
    SET_VAR(phi_rl_deg_vals);
    SET_VAR(rho_dot_kms_vals);
    SET_VAR(t_oet_spm_vals);
    SET_VAR(t_ret_spm_vals);
    SET_VAR(t_set_spm_vals);
    SET_VAR(rx_km_vals);
    SET_VAR(ry_km_vals);
    SET_VAR(rz_km_vals);
    SET_VAR(raw_tau_threshold_vals);
    SET_VAR(rho_corr_pole_km_vals);
    SET_VAR(rho_corr_timing_km_vals);
    SET_VAR(tau_power);
    SET_VAR(tau_phase);
    SET_VAR(tau_vals);
}
/*  End of crssringoccs_C_CSV_to_Py_CSV.                                      */
