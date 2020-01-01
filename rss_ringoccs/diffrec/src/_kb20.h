/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS__KB20_H
#define RSS_RINGOCCS__KB20_H

#include "special_functions.h"

static PyObject *kb20(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *x;
    double dx;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!d", &PyArray_Type, &x, &dx))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.kb20\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\tx:     Numpy Array of real numbers (Floats)\n"
            "\r\tdx:    Positive real number (Float)\n"
            "\rNotes:\n"
            "\r\tx must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(x);
    dim     = PyArray_DIMS(x)[0];
    data    = PyArray_DATA(x);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(x) != 1)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.kb20\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (dim == 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.kb20\n\n"
            "\rInput numpy array is empty.\n"
        );
        return NULL;
    }

    /*  Check that dx is positive.                                            */
    if (dx <= 0)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.kb20\n\n"
            "\rdx must be a positive number.\n"
        );
        return NULL;
    }

    if (typenum == NPY_FLOAT)
    {
        float *y;
        y = (float *)malloc(dim*sizeof(float));
        __get_one_real_from_two_real(((float *)data), dx, y, dim,
                                     Kaiser_Bessel_2_0_Float);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE)
    {
        long double *y;
        y = (long double *)malloc(dim*sizeof(long double));
        __get_one_real_from_two_real(((long double *)data), dx, y, dim,
                                     Kaiser_Bessel_2_0_Long_Double);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_LONGDOUBLE, (void *)y);
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);
    }
    else
    {
        double *y;
        y = (double *)malloc(dim*sizeof(double));

        if (typenum == NPY_DOUBLE)
        {
            __get_one_real_from_two_real(((double *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Double);
        }
        else if (typenum == NPY_BYTE)
        {
            __get_one_real_from_two_real(((char *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Char);
        }
        else if (typenum == NPY_UBYTE)
        {
            __get_one_real_from_two_real(((unsigned char *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_UChar);
        }
        else if (typenum == NPY_SHORT)
        {
            __get_one_real_from_two_real(((short *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Short);
        }
        else if (typenum == NPY_USHORT)
        {
            __get_one_real_from_two_real(((unsigned short *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_UShort);
        }
        else if (typenum == NPY_INT)
        {
            __get_one_real_from_two_real(((int *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Int);
        }
        else if (typenum == NPY_UINT)
        {
            __get_one_real_from_two_real(((unsigned int *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_UInt);
        }
        else if (typenum == NPY_LONG)
        {
            __get_one_real_from_two_real(((long *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Long);
        }
        else if (typenum == NPY_ULONG)
        {
            __get_one_real_from_two_real(((unsigned long *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_ULong);
        }
        else if (typenum == NPY_LONGLONG)
        {
            __get_one_real_from_two_real(((long long *)data), dx, y, dim,
                                         Kaiser_Bessel_2_0_Long_Long);
        }
        else if (typenum == NPY_ULONG)
        {
            __get_one_real_from_two_real(((unsigned long long *)data), dx, y,
                                         dim, Kaiser_Bessel_2_0_ULong_Long);
        }
        else
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.kb20\n\n"
                "\rInvalid data type for input numpy array. Input should be\n"
                "\ra one dimensional numpy array of real numbers (float).\n"
            );
            return NULL;
        }

        output = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)y);
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);
    }

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

#endif