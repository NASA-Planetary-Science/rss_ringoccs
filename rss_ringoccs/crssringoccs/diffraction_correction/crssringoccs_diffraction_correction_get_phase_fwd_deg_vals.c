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
 *  Author:     Ryan Maguire                                                  *
 *  Date:       March 10, 2026                                                *
 ******************************************************************************/

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  TMPL_MALLOC macro found here, providing C vs. C++ compatibility.          */
#include <libtmpl/include/compat/tmpl_malloc.h>

/*  Complex numbers and functions given here.                                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Radians to degrees conversion found here.                                 */
#include <libtmpl/include/constants/tmpl_math_constants.h>

/*  NULL macro found here.                                                    */
#include <stddef.h>

PyObject *
crssringoccs_DiffractionCorrection_Get_Phase_Fwd_Deg_Vals(PyObject *op,
                                                          void *closure)
{
    /*  Variable for the diffraction correction class instance.               */
    crssringoccs_PyDiffrecObj *self;

    /*  Buffer for the forward model phase data.                              */
    double *model_phase = NULL;

    /*  Variable for indexing over the array.                                 */
    size_t n;

    /*  If the input is NULL, there's nothing to be done.                     */
    if (!op)
        return NULL;

    /*  Get a pointer to the actual DiffractionCorrection instance.           */
    self = (crssringoccs_PyDiffrecObj *)op;

    /*  If the Tau variable has not been initialized, there is no data to     *
     *  process. Return None in this case.                                    */
    if (!self->tau)
    {
        MAKE_NONE(self->phase_fwd_deg_vals);
        Py_INCREF(self->phase_fwd_deg_vals);
        return self->phase_fwd_deg_vals;
    }

    /*  The real-valued phase_fwd_deg_vals array is computed from the complex *
     *  T_fwd array. If T_fwd is NULL, set phase_fwd_deg_vals to None.        */
    if (!self->tau->T_fwd)
    {
        MAKE_NONE(self->phase_fwd_deg_vals);
        Py_INCREF(self->phase_fwd_deg_vals);
        return self->phase_fwd_deg_vals;
    }

    /*  If we have already computed phase_fwd_deg_vals, increment the         *
     *  reference counter and return this object to the caller.               */
    if (self->phase_fwd_deg_vals)
    {
        Py_INCREF(self->phase_fwd_deg_vals);
        return self->phase_fwd_deg_vals;
    }

    /*  Otherwise, compute the model phase from the complex array.            */
    model_phase = TMPL_MALLOC(double, self->tau->n_used);

    /*  Check to make sure malloc did not fail.                               */
    if (!model_phase)
        return PyErr_NoMemory();

    /*  Otherwise, we have memory for the data. Compute the phase.            */
    for (n = 0; n < self->tau->n_used; ++n)
    {
        /*  The tau variables are shifted over to the start of the requested  *
         *  range for processing. Shift the index for T_fwd and get the data. */
        const size_t index = self->tau->start + n;
        const tmpl_ComplexDouble transmittance = self->tau->T_fwd[index];

        /*  We have T_hat = sqrt(power) * exp(i phase). From this, the phase  *
         *  is given by the argument of T_fwd. Compute this.                  */
        const double phase_rad = tmpl_CDouble_Argument(transmittance);

        /*  Convert from radians to degrees.                                  */
        model_phase[n] = phase_rad * tmpl_double_rad_to_deg;
    }

    /*  Wrap the data in a numpy array for ease of use.                       */
    crssringoccs_Create_Real_Numpy_Array(
        &self->phase_fwd_deg_vals,      /*  The Python object.                */
        model_phase,                    /*  The real-valued C data.           */
        crssringoccs_Capsule_Cleanup,   /*  Cleanup function for freeing data.*/
        self->tau->n_used               /*  Number of points in the array.    */
    );

    /*  We are returning a reference to this new array. Increment the counter.*/
    Py_INCREF(self->phase_fwd_deg_vals);
    return self->phase_fwd_deg_vals;
}
/*  End of crssringoccs_DiffractionCorrection_Get_Phase_Fwd_Deg_Vals.         */
