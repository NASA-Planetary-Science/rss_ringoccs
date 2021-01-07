#ifndef __RSS_RINGOCCS_AUXILIARY_H__
#define __RSS_RINGOCCS_AUXILIARY_H__

#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include "rss_ringoccs_py_api.h"

/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level.                  */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

#define _convert_void_to_double(in, out, n, dim, type)                         \
    out = malloc(sizeof(*out) * dim);                                          \
    for (n=0; n<dim; ++n)                                                      \
        out[n] = (double)((type *)in)[n];

static void *_get_array_from_one_long(void *in, long dim, long (*f)(long))
{
    long n, x;
    void *out;
    long *out_data;
    out_data = malloc(sizeof(*out_data)*dim);

    for (n=0; n<dim; ++n)
    {
        x = ((long *)in)[n];
        out_data[n] = f(x);
    }

    out = out_data;
    return out;
}

static void *_get_array_from_one_float(void *in, long dim, float (*f)(float))
{
    long n;
    void *out;
    float *out_data;
    float x;
    out_data = malloc(sizeof(*out_data) * dim);

    for (n=0; n<dim; ++n)
    {
        x = ((float *)in)[n];
        out_data[n] = f(x);
    }

    out = out_data;
    return out;
}

static void *_get_array_from_one_double(void *in, long dim, double (*f)(double))
{
    long n;
    void *out;
    double *out_data;
    out_data = malloc(sizeof(*out_data)*dim);

    for (n=0; n<dim; ++n)
        out_data[n] = f(((double *)in)[n]);

    out = out_data;
    return out;
}

static void *_get_array_from_one_ldouble(void *in, long dim,
                                         long double (*f)(long double))
{
    long n;
    void *out;
    long double *out_data;
    out_data = malloc(sizeof(*out_data)*dim);

    for (n=0; n<dim; ++n)
        out_data[n] = f(((long double *)in)[n]);

    out = out_data;
    return out;
}

static void *_get_complex_from_double(double *in, long dim,
                                      rssringoccs_ComplexDouble (*f)(double))
{
    long n;
    void *out;
    rssringoccs_ComplexDouble *out_data;
    out_data = malloc(sizeof(*out_data)*dim);

    for (n=0; n<dim; ++n)
        out_data[n] = f(in[n]);

    out = out_data;
    return out;
}

static void *
_get_complex_from_complex(void *in, long dim,
                          rssringoccs_ComplexDouble
                          (*f)(rssringoccs_ComplexDouble))
{
    long n;
    void *out;
    rssringoccs_ComplexDouble *out_data;
    rssringoccs_ComplexDouble z;
    out_data = malloc(sizeof(*out_data) * dim);

    for (n=0; n<dim; ++n)
    {
        z = ((rssringoccs_ComplexDouble *)in)[n];
        out_data[n] = f(z);
    }

    out = out_data;
    return out;
}

static void *_get_void_from_double(double *in, long dim, double (*f)(double))
{
    long n;
    void *out;
    double *out_data;
    out_data = malloc(sizeof(*out_data) * dim);

    for (n=0; n<dim; ++n)
        out_data[n] = f(in[n]);

    out = out_data;
    return out;
}

static void *_get_cvoid_from_double(double *in, long dim,
                                    rssringoccs_ComplexDouble (*f)(double))
{
    long n;
    void *out;
    rssringoccs_ComplexDouble *out_data;
    out_data = malloc(sizeof(*out_data) * dim);

    for (n=0; n<dim; ++n)
        out_data[n] = f(in[n]);

    out = out_data;
    return out;
}

static PyObject *
rssringoccs_Get_Out_From_Int(PyObject *x,
                             rssringoccs_Generic_Function_Obj *c_func)
{
    double x_val, y_val;
    long x_int, y_int;
    rssringoccs_ComplexDouble y_complex;

    if (c_func->long_func != NULL)
    {
        x_int = PyLong_AsLong(x);
        y_int = c_func->long_func(x_int);
        return PyLong_FromLong(y_int);
    }
    else if (c_func->double_func != NULL)
    {
        x_val = PyLong_AsDouble(x);
        y_val = c_func->double_func(x_val);
        return PyFloat_FromDouble(y_val);
    }
    else if (c_func->cdouble_from_real_func != NULL)
    {
        x_val     = PyLong_AsDouble(x);
        y_complex = c_func->cdouble_from_real_func(x_val);
        x_val = rssringoccs_CDouble_Real_Part(y_complex);
        y_val = rssringoccs_CDouble_Imag_Part(y_complex);
        return PyComplex_FromDoubles(x_val, y_val);
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\t%s\n\n"
                    "\rInteger input provided but this function does not\n"
                    "\rexcept real valued (float or int) arguments.",
                    c_func->func_name);
        return NULL;
    }
}

static PyObject *
rssringoccs_Get_Out_From_Float(PyObject *x,
                               rssringoccs_Generic_Function_Obj *c_func)
{
    double x_val, y_val;
    rssringoccs_ComplexDouble y_complex;

    if (c_func->double_func != NULL)
    {
        x_val = PyFloat_AsDouble(x);
        y_val = c_func->double_func(x_val);
        return PyFloat_FromDouble(y_val);
    }
    else if (c_func->cdouble_from_real_func != NULL)
    {
        x_val     = PyFloat_AsDouble(x);
        y_complex = c_func->cdouble_from_real_func(x_val);
        x_val = rssringoccs_CDouble_Real_Part(y_complex);
        y_val = rssringoccs_CDouble_Imag_Part(y_complex);
        return PyComplex_FromDoubles(x_val, y_val);
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\t%s\n\n"
                    "\rFloat input provided but this function does not\n"
                    "\rexcept real valued (float or int) arguments.",
                    c_func->func_name);
        return NULL;
    }
}

static PyObject *
rssringoccs_Get_Out_From_Complex(PyObject *x,
                                 rssringoccs_Generic_Function_Obj *c_func)
{
    double x_val, y_val;
    rssringoccs_ComplexDouble x_complex, y_complex;

    if (c_func->cdouble_from_complex_func != NULL)
    {
        x_val = PyComplex_RealAsDouble(x);
        y_val = PyComplex_ImagAsDouble(x);
        x_complex = rssringoccs_CDouble_Rect(x_val, y_val);
        y_complex = c_func->cdouble_from_complex_func(x_complex);
        x_val = rssringoccs_CDouble_Real_Part(y_complex);
        y_val = rssringoccs_CDouble_Imag_Part(y_complex);
        return PyComplex_FromDoubles(x_val, y_val);
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\t%s\n\n"
                    "\rFloat input provided but this function does not\n"
                    "\rexcept real valued (float or int) arguments.",
                    c_func->func_name);
        return NULL;
    }
}

static PyObject *
rssringoccs_Get_Py_Func_From_C(PyObject *self, PyObject *args,
                               rssringoccs_Generic_Function_Obj *c_func)
{
    /*  Declare necessary variables.                                          */
    PyObject *capsule, *output, *x, *nth_item;
    PyArrayObject *x_as_arr;
    long n, dim;
    double *new_x;
    rssringoccs_Bool error_occured = rssringoccs_False;
    void *data, *out;
    int typenum;
    out = NULL;

    if (c_func->func_name == NULL)
    {
        PyErr_Format(PyExc_RuntimeError,
                     "\n\rError Encountered: rss_ringoccs\n"
                     "\r\trssringoccs_Get_Py_Func_From_C\n\n"
                     "\rInput rssringoccs_Generic_Function_Obj does not\n"
                     "\rcontain a function name.\n\n");
        return NULL;
    }

    /*  Parse the data from Python and try to convert it to a usable format.  */
    if (!PyArg_ParseTuple(args, "O", &x))
    {
        PyErr_Format(PyExc_RuntimeError,
                     "\n\rError Encountered: rss_ringoccs\n"
                     "\r\t%s\n\n"
                     "\rCould not parse inputs.\n",
                     c_func->func_name);
        return NULL;
    }
    if (PyLong_Check(x))
        return rssringoccs_Get_Out_From_Int(x, c_func);
    else if (PyFloat_Check(x))
        return rssringoccs_Get_Out_From_Float(x, c_func);
    else if (PyComplex_Check(x))
        return rssringoccs_Get_Out_From_Complex(x, c_func);
    else if (PyList_Check(x))
    {
        dim    = PyList_Size(x);
        output = PyList_New(dim);

        for (n=0; n<dim; ++n)
        {
            nth_item = PyList_GET_ITEM(x, n);
            if (PyLong_Check(nth_item))
                PyList_SET_ITEM(output, n,
                                rssringoccs_Get_Out_From_Int(x, c_func));
            else if (PyFloat_Check(nth_item))
                PyList_SET_ITEM(output, n,
                                rssringoccs_Get_Out_From_Float(x, c_func));
            else if (PyComplex_Check(nth_item))
                PyList_SET_ITEM(output, n,
                                rssringoccs_Get_Out_From_Complex(x, c_func));
            else
            {
                PyErr_Format(PyExc_RuntimeError,
                             "\n\rError Encountered: rss_ringoccs\n"
                             "\r\t%s\n\n"
                             "\rInput list must contain real or\n"
                             "\rcomplex numbers only.\n",
                             c_func->func_name);
                return NULL;
            }
        }
        return output;
    }
    else if (!(PyArray_Check(x)))
    {
        PyErr_Format(PyExc_RuntimeError,
             "\n\rError Encountered: rss_ringoccs\n"
             "\r\t%s\n\n"
             "\rCould not parse inputs.\n",
             c_func->func_name);
        return NULL;
    }

    /*  If you get here, then the input is a numpy array. Grab some useful    *
     *  information about the data using numpy's API functions.               */
    x_as_arr = (PyArrayObject *)x;
    typenum = PyArray_TYPE(x_as_arr);
    dim     = (unsigned long)PyArray_DIM(x_as_arr, 0);
    data    = PyArray_DATA(x_as_arr);

    /*  Check the inputs to make sure they're valid.                         */
    if (PyArray_NDIM(x_as_arr) != 1)
    {
        PyErr_Format(PyExc_RuntimeError,
                     "\n\rError Encountered: rss_ringoccs\n"
                     "\r\tspecial_functions.%s\n"
                     "\n\rInput is not 1-dimensional.\n",
                     c_func->func_name);
        return NULL;
    }
    else if (dim == 0)
    {
        PyErr_Format(PyExc_RuntimeError,
                     "\n\rError Encountered: rss_ringoccs\n"
                     "\r\tspecial_functions.%s"
                     "\n\n\rInput numpy array is empty.\n",
                     c_func->func_name);
        return NULL;
    }

    switch (typenum)
    {
        case NPY_FLOAT:
            if (c_func->float_func != NULL)
                 out = _get_array_from_one_float(data, dim, c_func->float_func);
            else if (c_func->double_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, float);
                out = _get_void_from_double(new_x, dim, c_func->double_func);
                typenum = NPY_DOUBLE;
            }
            else if (c_func->cdouble_from_real_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, float);
                out = _get_cvoid_from_double(new_x, dim,
                                             c_func->cdouble_from_real_func);
                typenum = NPY_CDOUBLE;
            }
            else
                error_occured = rssringoccs_True;

            break;
        case NPY_DOUBLE:
            if (c_func->double_func != NULL)
                out = _get_array_from_one_double(data, dim,
                                                 c_func->double_func);
            else if (c_func->cdouble_from_real_func != NULL)
                out = _get_complex_from_double(data, dim,
                                               c_func->cdouble_from_real_func);
            else
                error_occured = rssringoccs_True;

            break;
        case NPY_LONGDOUBLE:
            if (c_func->ldouble_func != NULL)
                out = _get_array_from_one_ldouble(data, dim,
                                                  c_func->ldouble_func);
            else if (c_func->double_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, long double);
                out = _get_void_from_double(new_x, dim, c_func->double_func);
                typenum = NPY_DOUBLE;
            }
            else if (c_func->cdouble_from_real_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, long double);
                out = _get_cvoid_from_double(new_x, dim,
                                             c_func->cdouble_from_real_func);
                typenum = NPY_CDOUBLE;
            }
            else
                error_occured = rssringoccs_True;

            break;
        case NPY_CDOUBLE:
            if (c_func->cdouble_from_complex_func != NULL)
                out = _get_complex_from_complex(
                    data, dim, c_func->cdouble_from_complex_func
                );

            break;
        case NPY_LONG:
            if (c_func->long_func != NULL)
                out = _get_array_from_one_long(data, dim, c_func->long_func);
            else if (c_func->double_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, long);
                out = _get_void_from_double(new_x, dim, c_func->double_func);
                typenum = NPY_DOUBLE;
            }
            else if (c_func->cdouble_from_real_func != NULL)
            {
                _convert_void_to_double(data, new_x, n, dim, long);
                out = _get_cvoid_from_double(new_x, dim,
                                             c_func->cdouble_from_real_func);
                typenum = NPY_CDOUBLE;
            }
            else
                error_occured = rssringoccs_True;

        default:
            error_occured = rssringoccs_True;

    }

    if (error_occured)
    {
        PyErr_Format(PyExc_RuntimeError,
             "\n\rError Encountered: rss_ringoccs\n"
             "\r\t%s\n\n"
             "\rCould not parse inputs.\n",
             c_func->func_name);
        return NULL;
    }
    else if (out == NULL)
    {
        PyErr_Format(PyExc_RuntimeError,
             "\n\rError Encountered: rss_ringoccs\n"
             "\r\t%s\n\n"
             "\rCould not parse inputs. out returned NULL.\n",
             c_func->func_name);
        return NULL;
    }

    output  = PyArray_SimpleNewFromData(1, &dim, typenum, out);
    capsule = PyCapsule_New(out, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.  */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                     */
    return Py_BuildValue("N", output);
}

#endif
