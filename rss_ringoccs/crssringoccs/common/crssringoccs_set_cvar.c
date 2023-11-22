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
 *                        Diffraction Correction Class                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Defines the DiffractionCorrection class for rss_ringoccs. This uses   *
 *      the C-Python API to build an extension module containing the class.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/
#include "../crssringoccs.h"

/*  Avoid warnings about deprecated Numpy API versions.                       */
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  Creates a numpy array from a double array.                                */
void crssringoccs_Set_CVar(PyObject **py_ptr,
                           tmpl_ComplexDouble *ptr,
                           size_t len)
{
    PyObject *arr, *tmp, *capsule;
    npy_intp pylength = (npy_intp)len;

    /*  Numpy's _import_array function must be called before using the API.   */
    if (PyArray_API == NULL)
    {
        /*  If the import fails we can't safe use the tools. Abort.           */
        if (_import_array() < 0)
        {
            PyErr_Print();
            PyErr_SetString(PyExc_ImportError,
                            "numpy.core.multiarray failed to import");
            return;
        }
    }

    /*  If the pointer has memory allocated to it, create a numpy array.      */
    if (ptr != NULL)
    {
        /*  Numpy API function for creating numpy arrays from existing data.  */
        arr = PyArray_SimpleNewFromData(1, &pylength, NPY_CDOUBLE, ptr);

        /*  Create a capsule for this pointer so it is free'd when the numpy  *
         *  array is destroyed. Avoids memory leaks for the end-user.         */
        capsule = PyCapsule_New(ptr, NULL, crssringoccs_Capsule_Cleanup);

        /*  Link the array to the capsule. "del arr" in Python now free's the *
         *  memory allocated for the C pointer.                               */
        if (PyArray_SetBaseObject((PyArrayObject *)arr, capsule) == -1)
        {
            Py_DECREF(arr);
            return;
        }

        *py_ptr = Py_BuildValue("N", arr);
    }

    /*  Otherwise set the variable to a "None" object.                        */
    else
    {
        tmp = *py_ptr;
        Py_INCREF(Py_None);
        *py_ptr = Py_None;
        Py_XDECREF(tmp);
    }
}
/*  End of crssringoccs_set_var.                                              */
