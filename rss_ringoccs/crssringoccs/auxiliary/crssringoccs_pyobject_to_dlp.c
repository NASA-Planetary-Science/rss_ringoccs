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
#include "../crssringoccs_numpy_api.h"

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  NULL macro is given here.                                                 */
#include <stddef.h>

#define EXTRACT_VAR(var) \
    dlp->var = crssringoccs_DLP_Extract_Data(dlp, object, #var)

rssringoccs_DLPObj *crssringoccs_PyObject_To_DLP(PyObject * const object)
{
    PyObject *tmp = NULL;
    PyArrayObject *arr = NULL;
    rssringoccs_DLPObj *dlp = NULL;

    if (!object)
        return NULL;

    dlp = malloc(sizeof(*dlp));

    if (!dlp)
        return dlp;

    rssringoccs_DLP_Init(dlp);

    if (!object)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\rInput PyObject is NULL.\n";

        return dlp;
    }

    /*  Next we're going to run error checks on the input numpy arrays which  *
     *  should be contained inside of the DLPInst object. We'll check that    *
     *  these attributes exist, that they are numpy arrays, are 1 dimensional,*
     *  and have the same number of elements as rho_km_vals. We'll also       *
     *  convert the arrays to double and retrieve a pointer to the data.      *
     *  First, we need to make sure rho_km_vals is a legal numpy array and    *
     *  extract the length of it. Check that rho_km_vals exists in DLPInst.   */
    if (!PyObject_HasAttrString(object, "rho_km_vals"))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\trho_km_vals\n\n";

        return dlp;
    }

    /*  If it exists, get a pointer to it.                                    */
    tmp = PyObject_GetAttrString(object, "rho_km_vals");

    /*  Now make sure rho_km_vals is a numpy array.                           */
    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\rrho_km_vals must be a numpy array.\n";

        return dlp;
    }

    /*  If rho_km_vals is a numpy array, try to convert it to double.         */
    arr = (PyArrayObject *)tmp;

    if (PyArray_TYPE(arr) != NPY_DOUBLE)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\r\trho_km_vals must be an array of floats.\n\n";

        Py_CLEAR(tmp);
        return NULL;
    }

    if (PyArray_NDIM(arr) != 1)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\r\trho_km_vals must be a one dimensional numpy array.\n\n";

        Py_CLEAR(tmp);
        return NULL;
    }

    if (!PyArray_ISCARRAY(arr))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_PyObject_To_DLP\n\n"
            "\r\trho_km_vals must be a contiguous (C-style) numpy array.\n\n";

        Py_CLEAR(tmp);
        return NULL;
    }

    /*  If every passed, set tau.rho_km_vals to point to the data inside arr. */
    dlp->rho_km_vals = (double *)PyArray_DATA(arr);
    dlp->arr_size = PyArray_DIMS(arr)[0];

    Py_CLEAR(tmp);

    EXTRACT_VAR(p_norm_vals);
    EXTRACT_VAR(phase_deg_vals);
    EXTRACT_VAR(phi_deg_vals);
    EXTRACT_VAR(phi_rl_deg_vals);
    EXTRACT_VAR(B_deg_vals);
    EXTRACT_VAR(D_km_vals);
    EXTRACT_VAR(f_sky_hz_vals);
    EXTRACT_VAR(rho_dot_kms_vals);
    EXTRACT_VAR(t_oet_spm_vals);
    EXTRACT_VAR(t_ret_spm_vals);
    EXTRACT_VAR(t_set_spm_vals);
    EXTRACT_VAR(rx_km_vals);
    EXTRACT_VAR(ry_km_vals);
    EXTRACT_VAR(rz_km_vals);
    EXTRACT_VAR(rho_corr_pole_km_vals);
    EXTRACT_VAR(rho_corr_timing_km_vals);
    EXTRACT_VAR(raw_tau_threshold_vals);
    dlp->dx_km = dlp->rho_km_vals[1] - dlp->rho_km_vals[0];

    return rssringoccs_DLP_New_Reference(dlp);
}
