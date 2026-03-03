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

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  NULL macro defined here.                                                  */
#include <stddef.h>

PyObject *
crssringoccs_DiffractionCorrection_New(PyTypeObject *type,
                                       PyObject *args,
                                       PyObject *kwds)
{
    crssringoccs_PyDiffrecObj *self;
    PyObject *obj = type->tp_alloc(type, 0);

    if (!obj)
        return NULL;

    self = (crssringoccs_PyDiffrecObj *)obj;

    /*  The DiffractionCorrection class contains a pointer to a C struct, the *
     *  C-equivalent of the DiffractionCorrection class. This is not directly *
     *  accessible to Python users, but all of the data can be accessed via   *
     *  the other attributes (below) in the DiffractionCorrection class. Each *
     *  function that uses the Tau object checks for NULL pointers before     *
     *  trying to access. Initialize tau to NULL to avoid future errors.      */
    self->tau = NULL;

    /*  Set all of the arrays to NULL. These are populated with data by the   *
     *  init function once a DLP object is passed to DiffractionCorrection.   */
    self->T_in = NULL;
    self->T_out = NULL;
    self->T_fwd = NULL;
    self->k_vals = NULL;
    self->B_deg_vals = NULL;
    self->D_km_vals = NULL;
    self->F_km_vals = NULL;
    self->phi_deg_vals = NULL;
    self->phi_rl_deg_vals = NULL;
    self->rho_corr_pole_km_vals = NULL;
    self->rho_corr_timing_km_vals = NULL;
    self->rho_dot_kms_vals = NULL;
    self->rho_km_vals = NULL;
    self->t_oet_spm_vals = NULL;
    self->t_ret_spm_vals = NULL;
    self->t_set_spm_vals = NULL;
    self->tau_threshold_vals = NULL;
    self->w_km_vals = NULL;
    self->rx_km_vals = NULL;
    self->ry_km_vals = NULL;
    self->rz_km_vals = NULL;

    /*  The DiffractionCorrection class has a few more attributes that are    *
     *  Python objects. These can be lists, dictionaries, or strings. Each of *
     *  these are set when init is called, initialize them to NULL.           */
    self->input_vars = NULL;
    self->input_kwds = NULL;
    self->rngreq = NULL;
    self->perturb = NULL;

    /*  The kbmd20 is a new window, a modifed Kaiser-Bessel with alpha set to *
     *  two pi. The modification ensures the window goes to zero at its edges *
     *  while evaluating to one at the center, unlike the actual              *
     *  Kaiser-Bessel which is discontinuous at the edge of the window. The   *
     *  two pi factor is mostly guess work since it accurately reproduces the *
     *  PDS results. The real window used for that data is not known too me.  *
     *  The actual code for the window functions is in special_functions/     */
    self->wtype = "kbmd20";

    /*  Fresnel 4 is a new option, not mentioned in any of the papers but     *
     *  documented in our accompanying PDF. It uses Legendre polynomials to   *
     *  approximate the Fresnel kernel. It essentially takes Fresnels         *
     *  quadratic method to the next step, a quartic, hence the name. It is   *
     *  extremely fast (all of Rev007 takes less than a second) and very      *
     *  accurate for all but the most extreme occultations (like Rev133).     */
    self->psitype = "fresnel4";

    /*  By default, forward computations are not run, and the run is silent.  */
    self->use_fwd = tmpl_False;
    self->verbose = tmpl_False;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allan deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might be close to zero.                       */
    self->bfac = tmpl_True;
    self->use_norm = tmpl_True;

    /*  The default Allan deviation is the one for Cassini.                   */
    self->sigma = 2.0E-13;

    /*  Default resolution is 1 kilometer, same as the PDS.                   */
    self->input_resolution_km = 1.0;

    /*  If resolution_factor was not set, set to 0.75. This value was         *
     *  specified by Essam Marouf as necessary to ensure the reconstruction   *
     *  matches the PDS results. No justification is known to me.             */
    self->resolution_factor = 0.75;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->eccentricity = 0.0;
    self->periapse = 0.0;

    /*  Support for writing the TAB files is not yet available at the C level.*
     *  There are Python functions that handle this for now, set this to NULL.*/
    self->outfiles = NULL;

    /*  __new__ functions return PyObject points. Cast and return.            */
    return (PyObject *)self;
}
/*  End of crssringoccs_DiffractionCorrection_New.                            */
