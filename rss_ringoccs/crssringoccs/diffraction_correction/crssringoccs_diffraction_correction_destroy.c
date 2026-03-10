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

/*  Deallocating function for the DiffractionCorrection class.                */
void
crssringoccs_DiffractionCorrection_Destroy(crssringoccs_PyDiffrecObj *self)
{
    Py_CLEAR(self->T_in);
    Py_CLEAR(self->T_out);
    Py_CLEAR(self->T_fwd);
    Py_CLEAR(self->k_vals);
    Py_CLEAR(self->B_deg_vals);
    Py_CLEAR(self->D_km_vals);
    Py_CLEAR(self->F_km_vals);
    Py_CLEAR(self->phi_deg_vals);
    Py_CLEAR(self->phi_rl_deg_vals);
    Py_CLEAR(self->rho_corr_pole_km_vals);
    Py_CLEAR(self->rho_corr_timing_km_vals);
    Py_CLEAR(self->rho_dot_kms_vals);
    Py_CLEAR(self->rho_km_vals);
    Py_CLEAR(self->t_oet_spm_vals);
    Py_CLEAR(self->t_ret_spm_vals);
    Py_CLEAR(self->t_set_spm_vals);
    Py_CLEAR(self->tau_threshold_vals);
    Py_CLEAR(self->w_km_vals);
    Py_CLEAR(self->rx_km_vals);
    Py_CLEAR(self->ry_km_vals);
    Py_CLEAR(self->rz_km_vals);
    Py_CLEAR(self->input_vars);
    Py_CLEAR(self->input_kwds);
    Py_CLEAR(self->rngreq);
    Py_CLEAR(self->perturb);
    Py_CLEAR(self->p_norm_vals);
    Py_CLEAR(self->power_vals);
    Py_CLEAR(self->p_fwd_vals);
    Py_CLEAR(self->phase_norm_deg_vals);
    Py_CLEAR(self->phase_deg_vals);
    Py_CLEAR(self->phase_fwd_deg_vals);
    Py_CLEAR(self->tau_norm_vals);
    Py_CLEAR(self->tau_vals);
    Py_CLEAR(self->tau_fwd_vals);

    rssringoccs_Tau_Destroy(&self->tau);

    Py_TYPE(self)->tp_free((PyObject *) self);
}
/*  End of crssringoccs_DiffractionCorrection_Destroy.                        */
