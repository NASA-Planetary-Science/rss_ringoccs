#ifndef RSS_RINGOCCS__WINDOW_NORM_H
#define RSS_RINGOCCS__WINDOW_NORM_H

#include "special_functions.h"

static PyObject *window_norm(PyObject *self, PyObject *args){
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);
    double dx, f_scale;

    if (PyLong_Check(tuple)){
        long ker;
        PyArg_ParseTuple(args, "ldd", &ker, &dx, &f_scale);
        return PyFloat_FromDouble(
            Window_Normalization_Long(&ker, 1, dx, f_scale)
        );
    }
    else if (PyFloat_Check(tuple)){
        double ker;
        PyArg_ParseTuple(args, "ddd", &ker, &dx, &f_scale);
        return PyFloat_FromDouble(
            Window_Normalization_Double(&ker, 1, dx, f_scale)
        );
    }
    else if (PyArg_ParseTuple(args, "O!dd", &PyArray_Type,
                              &arr, &dx, &f_scale)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError,
                         "rss_ringoccs.diffrec.math_functions.min\n"
                         "\rInput must be a one-dimensional array.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (typenum == NPY_CFLOAT){
            return PyFloat_FromDouble(
                Window_Normalization_Complex_Float(
                    (complex float *)data, dim, dx, f_scale
                )
            );
        }
        else if (typenum == NPY_CDOUBLE){
            return PyFloat_FromDouble(
                Window_Normalization_Complex_Double(
                    (complex double *)data, dim, dx, f_scale
                )
            );
        }
        else if (typenum == NPY_CLONGDOUBLE){
            return PyFloat_FromDouble(
                Window_Normalization_Complex_Long_Double(
                    (complex long double *)data, dim, dx, f_scale
                )
            );
        }
        else if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                Window_Normalization_Float((float *)data, dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                Window_Normalization_Double((double *)data, dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                Window_Normalization_Long_Double((long double *)data,
                                                 dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_INT){
            return PyLong_FromLong(
                Window_Normalization_Int((int *)data, dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(
                Window_Normalization_Long((long *)data, dim, dx, f_scale)
            );
        }
        else {
            PyErr_Format(PyExc_TypeError,
                        "rss_ringoccs.diffrec.math_functions.min\n"
                        "\rInput should be a numpy array of real numbers"
                        "or a floating point/integer value.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "rss_ringoccs.diffrec.math_functions.min\n"
                     "\rInput should be a numpy array of numbers,"
                     "or a floating point/integer value.");
        return NULL;
    }  
}

#endif