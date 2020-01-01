#ifndef RSS_RINGOCCS__FRESNEL_SIN_H
#define RSS_RINGOCCS__FRESNEL_SIN_H

#include "special_functions.h"

static PyObject *fresnel_sin(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *x;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &x)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_sin\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\tx: Numpy Array of real numbers (Floats)\n"
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
    if (PyArray_NDIM(x) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_sin\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_sin\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    if (typenum == NPY_FLOAT){
        float *y;
        y = (float *)malloc(dim*sizeof(float));
        __get_one_real_from_one_real(((float *)data), y, dim,
                                     Fresnel_Sine_Taylor_to_Asymptotic_Float);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE){
        long double *y;
        y = (long double *)malloc(dim*sizeof(long double));
        __get_one_real_from_one_real(
            ((long double *)data), y, dim,
            Fresnel_Sine_Taylor_to_Asymptotic_Long_Double
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_LONGDOUBLE, (void *)y);
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);
    }
    else {
        double *y;
        y = (double *)malloc(dim*sizeof(double));

        if (typenum == NPY_DOUBLE)
        {
            __get_one_real_from_one_real(
                ((double *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_Double
            );

        }
        else if (typenum == NPY_BYTE)
        {
            __get_one_real_from_one_real(
                ((char *)data), y, dim, Fresnel_Sine_Taylor_to_Asymptotic_Char
            );
        }
        else if (typenum == NPY_UBYTE)
        {
            __get_one_real_from_one_real(
                ((unsigned char *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_UChar
            );
        }
        else if (typenum == NPY_SHORT)
        {
            __get_one_real_from_one_real(
                ((short *)data), y, dim, Fresnel_Sine_Taylor_to_Asymptotic_Short
            );
        }
        else if (typenum == NPY_USHORT)
        {
            __get_one_real_from_one_real(
                ((unsigned short *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_UShort
            );
        }
        else if (typenum == NPY_INT)
        {
            __get_one_real_from_one_real(
                ((int *)data), y, dim, Fresnel_Sine_Taylor_to_Asymptotic_Int
            );
        }
        else if (typenum == NPY_UINT)
        {
            __get_one_real_from_one_real(
                ((unsigned int *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_UInt
            );
        }
        else if (typenum == NPY_LONG)
        {
            __get_one_real_from_one_real(
                ((long *)data), y, dim, Fresnel_Sine_Taylor_to_Asymptotic_Long
            );
        }
        else if (typenum == NPY_ULONG)
        {
            __get_one_real_from_one_real(
                ((unsigned long *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_ULong
            );
        }
        else if (typenum == NPY_LONGLONG)
        {
            __get_one_real_from_one_real(
                ((long long *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_Long_Long
            );
        }
        else if (typenum == NPY_ULONG)
        {
            __get_one_real_from_one_real(
                ((unsigned long long *)data), y, dim,
                Fresnel_Sine_Taylor_to_Asymptotic_Long_Long
            );
        }
        else
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.fresnel_sin\n\n"
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