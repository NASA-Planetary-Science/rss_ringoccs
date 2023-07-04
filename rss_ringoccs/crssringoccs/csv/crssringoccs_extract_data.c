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

/*  malloc and NULL are defined here.                                         */
#include <stdlib.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

double *
extract_data(rssringoccs_DLPObj *dlp, PyObject *py_dlp, const char *var_name)
{
    PyObject *tmp;
    PyObject *arr;
    unsigned long len;

    if (dlp == NULL)
        return NULL;

    if (dlp->error_occurred)
        return NULL;

    if (py_dlp == NULL)
        return NULL;

    if (PyArray_API == NULL)
    {
        if (_import_array() < 0)
        {
            PyErr_Print();
            PyErr_SetString(PyExc_ImportError,
                            "numpy.core.multiarray failed to import");
            return NULL;
        }
    }

    if (!PyObject_HasAttrString(py_dlp, var_name))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\t%s\n\n",
            var_name
        );
        return NULL;
    }
    else
        tmp = PyObject_GetAttrString(py_dlp, var_name);

    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    len = (unsigned long)PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error.  */
    if (!arr)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a one-dimensional numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  arr should have the same number of elements as rho_km_vals.           */
    else if (len != dlp->arr_size)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s and rho_km_vals have a different number of elements.\n",
            var_name
        );
        return NULL;
    }

    /*  If every passed, set ptr to point to the data inside the array arr.   */
    return (double *)PyArray_DATA((PyArrayObject *)arr);
}
