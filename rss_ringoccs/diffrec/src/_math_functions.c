#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

#include "__max.h"
#include "__min.h"

static PyObject *max(PyObject *self, PyObject *args)
{
    #define FNAME "rss_ringoccs.diffrec.math_functions.max\n"
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long max;
        PyArg_ParseTuple(args, "l", &max);
        return PyLong_FromLong(max);
    }
    else if (PyFloat_Check(tuple)){
        double max;
        PyArg_ParseTuple(args, "d", &max);
        return PyFloat_FromDouble(max);
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput must be a one-dimensional array.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0){
            PyErr_Format(PyExc_TypeError, FNAME "\rInput array is empty.");
            return NULL;
        }

        if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(Max_Double((double *)data, dim));
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(Max_Long((long *)data, dim));
        }
        else {
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput should be a numpy array of numbers.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError, FNAME
                     "\rInput should be a numpy array of numbers.");
        return NULL;
    }
}

static PyObject *min(PyObject *self, PyObject *args)
{
    #define FNAME "rss_ringoccs.diffrec.math_functions.max\n"
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long min;
        PyArg_ParseTuple(args, "l", &min);
        return PyLong_FromLong(min);
    }
    else if (PyFloat_Check(tuple)){
        double min;
        PyArg_ParseTuple(args, "d", &min);
        return PyFloat_FromDouble(min);
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput must be a one-dimensional array.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(Min_Double((double *)data, dim));
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(Min_Long((long *)data, dim));
        }
        else {
            PyErr_Format(PyExc_TypeError, FNAME
                        "\rInput should be a numpy array of numbers,"
                        "or a floating point/integer value.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError, FNAME
                     "\rInput should be a numpy array of numbers,"
                     "or a floating point/integer value.");
        return NULL;
    }
}

static PyMethodDef DiffMethods[] =
{
    {"max", max, METH_VARARGS, "Compute the maximum of a numpy array."},
    {"min", min, METH_VARARGS, "Compute the minimum of a numpy array."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cModPyDem =
    {PyModuleDef_HEAD_INIT, "_math_functions", "", -1, DiffMethods};

PyMODINIT_FUNC PyInit__math_functions(void)
{
    import_array();
    return PyModule_Create(&cModPyDem);
}