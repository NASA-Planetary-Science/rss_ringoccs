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

/*  NULL and size_t given here.                                               */
#include <stddef.h>

/*  Buffer for the error message that is stored in the DLP object should an   *
 *  error occur. Note that because we are using a global buffer, this routine *
 *  is not reentrent. Fortunately, the functions that call this routine are   *
 *  run in a single-threaded environment, so this causes no real issue.       */
static char rssringoccs_extract_data_error_message[1024];

double *
crssringoccs_DLP_Extract_Data(rssringoccs_DLPObj * const dlp,
                              PyObject * const object,
                              const char * const var_name)
{
    PyObject *tmp = NULL;
    PyArrayObject *arr = NULL;
    double *data = NULL;
    size_t len;

    if (!dlp)
        return NULL;

    if (dlp->error_occurred)
        return NULL;

    if (!object)
        return NULL;

    if (!PyObject_HasAttrString(object, var_name))
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\t%s\n\n",
            var_name
        );

        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    tmp = PyObject_GetAttrString(object, var_name);

    if (!tmp)
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\rPyObject_GetAttrString returned NULL for %s\n\n",
            var_name
        );

        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    if (!PyArray_Check(tmp))
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );

        Py_CLEAR(tmp);
        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    arr = (PyArrayObject *)tmp;

    if (PyArray_TYPE(arr) != NPY_DOUBLE)
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\r\t%s must be an array of floats.\n\n",
            var_name
        );

        Py_CLEAR(tmp);
        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    if (PyArray_NDIM(arr) != 1)
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\r\t%s must be a one dimensional numpy array.\n\n",
            var_name
        );

        Py_CLEAR(tmp);
        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    if (!PyArray_ISCARRAY(arr))
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\r\t%s must be a contiguous (C-style) numpy array.\n\n",
            var_name
        );

        Py_CLEAR(tmp);
        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    len = (size_t)PyArray_DIMS(arr)[0];

    /*  arr should have the same number of elements as rho_km_vals.           */
    if (len != dlp->arr_size)
    {
        sprintf(
            rssringoccs_extract_data_error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DLP_Extract_Data\n\n"
            "\r%s and rho_km_vals have a different number of elements.\n\n",
            var_name
        );

        dlp->error_occurred = tmpl_True;
        dlp->error_message = rssringoccs_extract_data_error_message;
        return NULL;
    }

    /*  If every passed, set ptr to point to the data inside the array arr.   */
    data = (double *)PyArray_DATA(arr);
    Py_CLEAR(tmp);
    return data;
}
