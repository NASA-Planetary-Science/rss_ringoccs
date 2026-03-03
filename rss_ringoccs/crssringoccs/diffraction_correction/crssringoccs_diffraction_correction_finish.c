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

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Macro for the crssringoccs_set_var function to shorten the syntax.        */
#define SET_REAL_VAR(a)                                                        \
crssringoccs_Create_Real_Numpy_Array(                                          \
    &self->a,                                                                  \
    self->tau->a + self->tau->start,                                           \
    NULL,                                                                      \
    self->tau->n_used                                                          \
)

#define SET_COMPLEX_VAR(a)                                                     \
crssringoccs_Create_Complex_Numpy_Array(                                       \
    &self->a,                                                                  \
    self->tau->a + self->tau->start,                                           \
    NULL,                                                                      \
    self->tau->n_used                                                          \
)

#define SET_DLP_VAR(a)                                                         \
    do {                                                                       \
        PyObject *tmp = PyObject_GetAttrString(dlp, #a);                       \
        self->a = PyObject_GetItem(tmp, slice);                                \
        Py_CLEAR(tmp);                                                         \
    } while (0)

/*  Macro for safely creating None objects.                                   */
#define MAKE_NONE(var)                                                         \
    do {                                                                       \
        PyObject *tmp = self->var;                                             \
        Py_INCREF(Py_None);                                                    \
        self->var = Py_None;                                                   \
        Py_XDECREF(tmp);                                                       \
    } while(0)

/*  Converts a C Tau struct to a Python Tau Object.                           */
void
crssringoccs_DiffractionCorrection_Finish(crssringoccs_PyDiffrecObj *self,
                                          PyObject *dlp)
{
    PyObject *slice = NULL;
    PyObject *start = NULL;
    PyObject *end = NULL;

    /*  If the C version of the object is NULL there is nothing to do.        */
    if (!self)
        return;

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    start = PyLong_FromSize_t(self->tau->start);
    end = PyLong_FromSize_t(self->tau->start + self->tau->n_used);
    slice = PySlice_New(start, end, NULL);

    if (!start)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Finish\n\n"
            "\rPyLong_FromSize_t returned NULL for start.\n\n";

        goto CLEANUP;
    }

    if (!end)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Finish\n\n"
            "\rPyLong_FromSize_t returned NULL for end.\n\n";

        goto CLEANUP;
    }

    if (!slice)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Finish\n\n"
            "\tPySlice_New returned NULL for slice.\n\n";

        goto CLEANUP;
    }

    /*  Set every variable in the Python object from the C Tau struct.        */
    SET_COMPLEX_VAR(T_in);
    SET_COMPLEX_VAR(T_out);

    SET_REAL_VAR(F_km_vals);
    SET_REAL_VAR(k_vals);
    SET_REAL_VAR(w_km_vals);
    MAKE_NONE(tau_threshold_vals);

    SET_DLP_VAR(rho_km_vals);
    SET_DLP_VAR(phi_deg_vals);
    SET_DLP_VAR(rho_dot_kms_vals);
    SET_DLP_VAR(B_deg_vals);
    SET_DLP_VAR(D_km_vals);
    SET_DLP_VAR(t_oet_spm_vals);
    SET_DLP_VAR(t_ret_spm_vals);
    SET_DLP_VAR(t_set_spm_vals);
    SET_DLP_VAR(rho_corr_pole_km_vals);
    SET_DLP_VAR(rho_corr_timing_km_vals);
    SET_DLP_VAR(phi_rl_deg_vals);
    SET_DLP_VAR(rx_km_vals);
    SET_DLP_VAR(ry_km_vals);
    SET_DLP_VAR(rz_km_vals);

    /*  If forward modeling was not performed, set these as None objects.     */
    if (!self->tau->T_fwd)
        MAKE_NONE(T_fwd);
    else
        SET_COMPLEX_VAR(T_fwd);

    CLEANUP:
        Py_CLEAR(slice);
        Py_CLEAR(start);
        Py_CLEAR(end);
}
/*  End of crssringoccs_C_Tau_To_Py_Tau.                                      */

/*  Undefine these in case someone wants to #include this file.               */
#undef SET_REAL_VAR
#undef SET_COMPLEX_VAR
#undef SET_DLP_VAR
#undef MAKE_NONE
