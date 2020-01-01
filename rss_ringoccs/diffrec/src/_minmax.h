/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS__MINMAX_H
#define RSS_RINGOCCS__MINMAX_H

#include "special_functions.h"

static PyObject *max(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long max;
        if (PyArg_ParseTuple(args, "l", &max)){
            return PyLong_FromLong(max);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                         "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double max;
        if(PyArg_ParseTuple(args, "d", &max)){
            return PyFloat_FromDouble(max);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                         "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1) goto FAIL;

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0) goto FAIL;

        switch(typenum)
        {
            case NPY_FLOAT:
                return PyFloat_FromDouble(Max_Float((float *)data, dim));
            case NPY_DOUBLE:
                return PyFloat_FromDouble(Max_Double((double *)data, dim));
            case NPY_LONGDOUBLE:
                return
                PyFloat_FromDouble(Max_Long_Double((long double *)data, dim));
            case NPY_SHORT:
                return PyLong_FromLong(Max_Short((short *)data, dim));
            case NPY_INT:
                return PyLong_FromLong(Max_Int((int *)data, dim));
            case NPY_LONG:
                return PyLong_FromLong(Max_Long((long *)data, dim));
            case NPY_LONGLONG:
                return PyLong_FromLong(Max_Long_Long((long long *)data, dim));
            default:
                goto FAIL;
        }
    }
    else goto FAIL;
    FAIL: {
        PyErr_Format(PyExc_TypeError,
                     "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                     "\r\t\tInput should be a one dimensional numpy array of\n"
                     "\r\t\treal numbers, or a float/int number.\n"
                     "\r\t\tExample:\n"
                     "\r\t\t\t>>> import numpy\n"
                     "\r\t\t\t>>> import _special_functions as sf\n"
                     "\r\t\t\t>>> x = numpy.random.rand(100)\n"
                     "\r\t\t\t>>> y = sf.max(x)\n\n"
                     "\r\t\tNOTE:\n"
                     "\r\t\t\tOnly one dimensional numpy arrays are allowed.\n"
                     "\r\t\t\tComplex numbers are not allowed. If the input\n"
                     "\r\t\t\tis a single floating point or integer number,\n"
                     "\r\t\t\tthe output will simply be that number.");
        return NULL;
    };
}

static PyObject *min(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long min;
        if (PyArg_ParseTuple(args, "l", &min)){
            return PyLong_FromLong(min);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double min;
        if(PyArg_ParseTuple(args, "d", &min)){
            return PyFloat_FromDouble(min);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput must be one dimensional.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0){
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput is zero dimensional.");
            return NULL;
        }

        if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                Min_Float((float *)data, dim)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                Min_Double((double *)data, dim)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                Min_Long_Double((long double *)data, dim)
            );
        }
        else if (typenum == NPY_SHORT){
            return PyLong_FromLong(
                Min_Short((short *)data, dim)
            );
        }
        else if (typenum == NPY_INT){
            return PyLong_FromLong(
                Min_Int((int *)data, dim)
            );
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(
                Min_Long((long *)data, dim)
            );
        }
        else if (typenum == NPY_LONGLONG){
            return PyLong_FromLong(
                Min_Long_Long((long long *)data, dim)
            );
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput should be a numpy array of numbers.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "rss_ringoccs.diffrec.math_functions.min\n"
                     "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                     "\r\t\tInput should be a numpy array of numbers.");
        return NULL;
    }
}

#endif