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
#include "../crssringoccs.h"

void ExtractCSVData_dealloc(PyCSVObj *self)
{
    Py_XDECREF(self->B_deg_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->f_sky_hz_vals);
    Py_XDECREF(self->p_norm_vals);
    Py_XDECREF(self->raw_tau_vals);
    Py_XDECREF(self->phase_deg_vals);
    Py_XDECREF(self->phi_deg_vals);
    Py_XDECREF(self->phi_rl_deg_vals);
    Py_XDECREF(self->raw_tau_threshold_vals);
    Py_XDECREF(self->rev_info);
    Py_XDECREF(self->rho_corr_pole_km_vals);
    Py_XDECREF(self->rho_corr_timing_km_vals);
    Py_XDECREF(self->rho_dot_kms_vals);
    Py_XDECREF(self->rho_km_vals);
    Py_XDECREF(self->t_oet_spm_vals);
    Py_XDECREF(self->t_ret_spm_vals);
    Py_XDECREF(self->t_set_spm_vals);
    Py_XDECREF(self->tau_vals);
    Py_XDECREF(self->history);
    Py_XDECREF(self->input_vars);
    Py_XDECREF(self->input_kwds);
    Py_XDECREF(self->rx_km_vals);
    Py_XDECREF(self->ry_km_vals);
    Py_XDECREF(self->rz_km_vals);
    Py_XDECREF(self->tau_phase);
    Py_XDECREF(self->tau_power);
    Py_XDECREF(self->tau_vals);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

