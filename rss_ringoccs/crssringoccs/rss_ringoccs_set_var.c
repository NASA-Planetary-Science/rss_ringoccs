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
#include <rss_ringoccs/include/crss_ringoccs.h>

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

void set_var(PyObject **py_ptr, double **ptr, unsigned long int len)
{
    PyObject *arr, *tmp, *capsule;
    long pylength = (long)len;

    if (PyArray_API == NULL)
    {
        if (_import_array() < 0)
        {
            PyErr_Print();
            PyErr_SetString(PyExc_ImportError,
                            "numpy.core.multiarray failed to import");
            return;
        }
    }

    if (*ptr != NULL)
    {
        arr = PyArray_SimpleNewFromData(1, &pylength, NPY_DOUBLE, *ptr);
        capsule = PyCapsule_New((void *) (*ptr), NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)arr, capsule);
        tmp = *py_ptr;
        Py_INCREF(arr);
        *py_ptr = arr;
        Py_XDECREF(tmp);
    }
    else
    {
        tmp = *py_ptr;
        Py_INCREF(Py_None);
        *py_ptr = Py_None;
        Py_XDECREF(tmp);
    }
}

