#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "../../include/Python.h"
#include "../../include/arrayobject.h"
#include <math.h>

static PyObject *MaxValue(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyArrayObject *arr2;
    int i;

    /*  parse single numpy array argument */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        PyErr_Format(PyExc_RuntimeError, "You done messed up.");
        return NULL;
    }

    // Useful metadata about the data
    npy_intp* dims = PyArray_DIMS(arr);
    npy_intp* strides = PyArray_STRIDES(arr);
    void *data0 = PyArray_DATA(arr);
    int typenum = PyArray_TYPE(arr);
    if( typenum == NPY_DOUBLE ){

    double max = *(double *)data0;

    for (i=0; i<n; ++i){
        data0 += strides[0];
        if (*(double *)data0 > max){
            max = *(double *)data0;
        }
    }

    /*
     *  Another way to do it:
     *      void* data0 = PyArray_DATA(arr);
     *      npy_intp *strides = PyArray_STRIDES (arr);
     *      double test, max = *(double *)&data0[0];
     *
     *      for (i=0; i<n; ++i){
     *          test = *(double *)&data0[i*strides[0]];
     *          if (test > max){
     *              max = test;
     *          }
     *      }
    */

    return Py_BuildValue("d", max);
}

/*  define functions in module */
static PyMethodDef DiffMethods[] =
{
    {
        "MaxValue",
        MaxValue,
        METH_VARARGS,
        "Evaluate the maximum of a numpy array."
    },
    {
        NULL,
        NULL,
        0,
        NULL
    }
};

static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "testmodule", /* name of module */
    "",
    -1,
    DiffMethods
};

/* module initialization */
PyMODINIT_FUNC PyInit_testmodule(void)
{
    /* IMPORTANT: this must be called */
    import_array();
    return PyModule_Create(&cModPyDem);
}
