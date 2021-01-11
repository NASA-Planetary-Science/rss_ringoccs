/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  The following are NON-STANDARD header files that MUST BE IN YOUR PATH.    *
 *  If you installed python using anaconda then Python.h should automatically *
 *  be included in your path. Also, if you are using the setup.py script      *
 *  provided then inclusion of these files should be done for you.            */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

#include <rss_ringoccs/include/rss_ringoccs_fft.h>

static PyObject *fft(PyObject *self, PyObject *args)
{
    /*  We'll need output and capsule for safely creating the output array    *
     *  and ensuring we don't have a memory leak. rho is the input variable.  */
    PyObject *output, *capsule, *arr_in, *arr;

    /*  Variable for the size of the input array or list.                     */
    long dim;
    int inverse;

    /*  Variables needed for the array data.                                  */
    rssringoccs_ComplexDouble *in, *out;

    /*  Try to parse the user input, returning error if this fails.           */
    if (!PyArg_ParseTuple(args, "Op", &arr_in, &inverse))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tspecial_functions.fft\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trx:      Numpy Array of real or complex numbers (Floats)\n\n"
            "\rKeywords:\n"
            "\r\tinverse: Compute inverse FFT.\n"
        );
        return NULL;
    }

    /*  If the user supplied a float or int, simply return the value.         */
    if (PyLong_Check(arr_in) || PyFloat_Check(arr_in))
        return arr_in;

    PyArray_Descr *arr_type = PyArray_DescrFromType(NPY_CDOUBLE);
    arr = PyArray_FromAny(arr_in, arr_type, 1, 1, NPY_ARRAY_BEHAVED, NULL);

    if (!arr)
    {
        PyErr_Format(PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tspecial_functions.fft\n\n"
                    "\rInvalid data type for one of the input arrays. Input\n"
                    "\rshoule be a 1-dimensional array of real numbers.\n");
        return NULL;
    }

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tspecial_functions.fft\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }

    /*  Get the size of the input numpy array.                                */
    dim = PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  Check that the array isn't empty. Raise error otherwise.              */
    if (dim == 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tspecial_functions.fft\n\n"
            "\rInput numpy array is empty.\n"
        );
        return NULL;
    }

    /*  Get a pointer to the actual data from the array. Allocate memory for  *
     *  the data of the output numpy array, which we'll call T_hat.           */
    in  = (rssringoccs_ComplexDouble *)PyArray_DATA((PyArrayObject *)arr);
    out = rssringoccs_Complex_FFT(in, dim, inverse);

    /*  Set the output and capsule, ensuring no memory leaks occur.           */
    output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)out);
    capsule = PyCapsule_New((void *)out, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}
