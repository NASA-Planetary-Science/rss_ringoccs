#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "../../include/Python.h"
#include "../../include/arrayobject.h"
#include <math.h>
#include <stdio.h>

static PyObject *max(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    long i, n, strides;

    if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        // Useful metadata about the data
        n = PyArray_DIMS(arr)[0];
        strides = PyArray_STRIDES(arr)[0];
        void *data0 = PyArray_DATA(arr);
        int typenum = PyArray_TYPE(arr);
        if(typenum == NPY_DOUBLE){
            double max = *(double *)data0;
            for (i=0; i<n; ++i){
                if (*(double *)data0 > max){
                    max = *(double *)data0;
                }
                data0 += strides;
            }
            return Py_BuildValue("d", max);
        }
        else if (typenum == NPY_LONG){
            long max = *(long *)data0;
            for (i=0; i<n; ++i){
                if (*(long *)data0 > max){
                    max = *(long *)data0;
                }
                data0 += strides;
            }
            return Py_BuildValue("l", max);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\rInput should be a numpy array of numbers (Floats/Ints)."
            );
            return NULL;
        }
    }
    else{
        PyErr_Format(
            PyExc_TypeError,
            "\rInput should be a numpy array of numbers (Floats/Ints)."
        );
        return NULL;
    }
}

/*  define functions in module */
static PyMethodDef DiffMethods[] =
{
    {"max", max, METH_VARARGS, "Evaluate the maximum of a numpy array."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "_math_functions", /* name of module */
    "",
    -1,
    DiffMethods
};

/* module initialization */
PyMODINIT_FUNC PyInit__math_functions(void)
{
    /* IMPORTANT: this must be called */
    import_array();
    return PyModule_Create(&cModPyDem);
}
