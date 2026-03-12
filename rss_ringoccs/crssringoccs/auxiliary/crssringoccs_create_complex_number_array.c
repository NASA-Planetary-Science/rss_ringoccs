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
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/
#include "../crssringoccs.h"
#include "../crssringoccs_numpy_api.h"

/*  Creates a numpy array from a double array.                                */
void
crssringoccs_Create_Complex_Numpy_Array(PyObject ** const py_ptr,
                                        tmpl_ComplexDouble * const ptr,
                                        void (*cleanup)(PyObject *),
                                        const size_t len)
{
    PyObject *arr = NULL;
    PyObject *capsule = NULL;
    npy_intp pylength = (npy_intp)len;

    /*  If the pointer has memory allocated to it, create a numpy array.      */
    if (ptr)
    {
        /*  Numpy API function for creating numpy arrays from existing data.  */
        arr = PyArray_SimpleNewFromData(1, &pylength, NPY_CDOUBLE, ptr);

        /*  Check if numpy failed to create a new array.                      */
        if (!arr)
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tcrssringoccs_Create_Complex_Numpy_Array\n\n"
                "\rPyArray_SimpleNewFromData returned NULL.\n\n"
            );

            *py_ptr = NULL;
            return;
        }

        /*  Check if we need to attach a memory cleanup function to the array.*/
        if (cleanup)
        {
            /*  Create a capsule for this pointer so it is free'd when the    *
             *  array is destroyed. Avoids memory leaks for the end-user.     */
            capsule = PyCapsule_New(ptr, NULL, cleanup);

            /*  Check if the new capsule is valid.                            */
            if (!capsule)
            {
                PyErr_Format(
                    PyExc_RuntimeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tcrssringoccs_Create_Complex_Numpy_Array\n\n"
                    "\rPyCapsule_New returned NULL.\n\n"
                );

                Py_CLEAR(arr);
                *py_ptr = NULL;
                return;
            }

            /*  Link the array to the capsule. "del arr" in Python now free's *
             *  the memory allocated for the C pointer.                       */
            if (PyArray_SetBaseObject((PyArrayObject *)arr, capsule) == -1)
            {
                PyErr_Format(
                    PyExc_RuntimeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tcrssringoccs_Create_Complex_Numpy_Array\n\n"
                    "\rPyArray_SetBaseObject failed to set capsule.\n\n"
                );

                Py_CLEAR(capsule);
                Py_CLEAR(arr);
                *py_ptr = NULL;
                return;
            }
        }

        *py_ptr = arr;
    }

    /*  Otherwise set the variable to a "None" object.                        */
    else
        MAKE_NONE(*py_ptr);
}
/*  End of crssringoccs_Create_Complex_Numpy_Array.                           */
