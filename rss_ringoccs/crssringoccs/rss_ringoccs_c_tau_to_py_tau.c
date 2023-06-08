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

void rssringoccs_C_Tau_to_Py_Tau(PyDiffrecObj *py_tau, rssringoccs_TAUObj *tau)
{
    PyObject *tmp;
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (py_tau == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput py_tau is NULL. Aborting.n"
        );
        return;
    }

    set_var(&py_tau->rho_km_vals, &tau->rho_km_vals, tau->arr_size);
    set_var(&py_tau->B_rad_vals, &tau->B_rad_vals, tau->arr_size);
    set_var(&py_tau->D_km_vals, &tau->D_km_vals, tau->arr_size);
    set_var(&py_tau->F_km_vals, &tau->F_km_vals, tau->arr_size);
    set_var(&py_tau->f_sky_hz_vals, &tau->f_sky_hz_vals, tau->arr_size);
    set_var(&py_tau->p_norm_vals, &tau->p_norm_vals, tau->arr_size);
    set_var(&py_tau->phase_rad_vals, &tau->phase_rad_vals, tau->arr_size);
    set_var(&py_tau->phase_vals, &tau->phase_vals, tau->arr_size);
    set_var(&py_tau->phi_rad_vals, &tau->phi_rad_vals, tau->arr_size);
    set_var(&py_tau->phi_rl_rad_vals, &tau->phi_rl_rad_vals, tau->arr_size);
    set_var(&py_tau->power_vals, &tau->power_vals, tau->arr_size);
    set_var(&py_tau->rho_dot_kms_vals, &tau->rho_dot_kms_vals, tau->arr_size);
    set_var(&py_tau->t_oet_spm_vals, &tau->t_oet_spm_vals, tau->arr_size);
    set_var(&py_tau->t_ret_spm_vals, &tau->t_ret_spm_vals, tau->arr_size);
    set_var(&py_tau->t_set_spm_vals, &tau->t_set_spm_vals, tau->arr_size);
    set_var(&py_tau->tau_vals, &tau->tau_vals, tau->arr_size);
    set_var(&py_tau->w_km_vals, &tau->w_km_vals, tau->arr_size);
    set_var(&py_tau->rx_km_vals, &tau->rx_km_vals, tau->arr_size);
    set_var(&py_tau->ry_km_vals, &tau->ry_km_vals, tau->arr_size);
    set_var(&py_tau->rz_km_vals, &tau->rz_km_vals, tau->arr_size);
    set_var(&py_tau->raw_tau_threshold_vals, &tau->raw_tau_threshold_vals, tau->arr_size);
    set_var(&py_tau->rho_corr_pole_km_vals, &tau->rho_corr_pole_km_vals, tau->arr_size);
    set_var(&py_tau->rho_corr_timing_km_vals, &tau->rho_corr_timing_km_vals, tau->arr_size);
    set_var(&py_tau->tau_threshold_vals, &tau->tau_threshold_vals, tau->arr_size);

    if (tau->T_fwd == NULL)
    {
        tmp = py_tau->p_norm_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->p_norm_fwd_vals = Py_None;
        Py_XDECREF(tmp);

        tmp = py_tau->phase_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->phase_fwd_vals = Py_None;
        Py_XDECREF(tmp);

        tmp = py_tau->tau_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->tau_fwd_vals = Py_None;
        Py_XDECREF(tmp);
    }
    else
    {
        set_var(&py_tau->p_norm_fwd_vals, &tau->p_norm_fwd_vals, tau->arr_size);
        set_var(&py_tau->phase_fwd_vals, &tau->phase_fwd_vals, tau->arr_size);
        set_var(&py_tau->tau_fwd_vals, &tau->tau_fwd_vals, tau->arr_size);
    }
}

/*  To edit:
    PyObject          *rev_info;
    PyObject          *history;
*/
