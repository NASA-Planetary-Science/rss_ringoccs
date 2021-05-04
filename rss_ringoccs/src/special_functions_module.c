/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  The following are NON-STANDARD header files that MUST BE IN YOUR PATH.    *
 *  If you installed python using anaconda then Python.h should automatically *
 *  be included in your path. Also, if you are using the setup.py script      *
 *  provided then inclusion of these files should be done for you.            */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  The following header files are NON-STANDARD, and are a part of the        *
 *  rss_ringoccs package. The setup scripts will add the correct CFLAGS so    *
 *  compiler should see these without user interaction.                       */

/*  This file contains a library, written with pure C99, of mathematical      *
 *  functions for user convenience, limit dependencies, and provide           *
 *  algorithms for the curious user which are accurate, but relatively        *
 *  simple. It contains Bessel, Fresnel, and Lambert functions, as well as    *
 *  some physics-based conversion functions (frequency-to-wavelength, etc.).  */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fft.h>

#include "auxiliary.h"

#define VarToString(Var) (#Var)

/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 *  This contains the Numpy-C and Python-C API parts that allow for the above *
 *  functions to be called in Python. Numpy arrays, as well as floating point *
 *  and integer valued arguments may then be passed into these functions for  *
 *  improvement in performance, as opposed to the routines written purely in  *
 *  Python. Successful compiling requires the Numpy and Python header files.  *
 *----------------------------------------------------------------------------*/

#ifdef RSS_RINGOCCS_GetArrFromTwoTypeFunc
#undef RSS_RINGOCCS_GetArrFromTwoTypeFunc
#endif

#define RSS_RINGOCCS_GetArrFromTwoTypeFunc(type, name)                         \
static void *_get_array_from_two_##name(void *in, type param1, long dim,       \
                                        type (*f)(type, type))                 \
{                                                                              \
    long n;                                                                    \
    void *out;                                                                 \
    type *out_data;                                                            \
                                                                               \
    out_data = malloc(sizeof(*out_data)*dim);                                  \
                                                                               \
    for (n=0; n<dim; ++n)                                                      \
        out_data[n] = f(((type *)in)[n], param1);                              \
                                                                               \
    out = out_data;                                                            \
    return out;                                                                \
}


RSS_RINGOCCS_GetArrFromTwoTypeFunc(float, float)
RSS_RINGOCCS_GetArrFromTwoTypeFunc(double, double)
RSS_RINGOCCS_GetArrFromTwoTypeFunc(long double, longdouble)

#undef RSS_RINGOCCS_GetArrFromTwoTypeFunc

#ifdef RSS_RINGOCCS_GetArrFromThreeTypeFunc
#undef RSS_RINGOCCS_GetArrFromThreeTypeFunc
#endif

#define RSS_RINGOCCS_GetArrFromThreeTypeFunc(type, name)                       \
static void *_get_array_from_three_##name(void *in, type param1,               \
                                          type param2, long dim,               \
                                          type (*f)(type, type, type))         \
{                                                                              \
    long n;                                                                    \
    void *out;                                                                 \
    type *out_data;                                                            \
                                                                               \
    out_data = malloc(sizeof(*out_data)*dim);                                  \
                                                                               \
    for (n=0; n<dim; ++n)                                                      \
        out_data[n] = f(((type *)in)[n], param1, param2);                      \
                                                                               \
    out = out_data;                                                            \
    return out;                                                                \
}

RSS_RINGOCCS_GetArrFromThreeTypeFunc(float, float)
RSS_RINGOCCS_GetArrFromThreeTypeFunc(double, double)
RSS_RINGOCCS_GetArrFromThreeTypeFunc(long double, longdouble)

#undef RSS_RINGOCCS_GetArrFromThreeTypeFunc

#define WindowFunctionForNumpy(FuncName, CName)                                \
static PyObject * FuncName(PyObject *self, PyObject *args)                     \
{                                                                              \
    PyObject *output, *capsule, *x, *nth_item;                                 \
    double W, x_val, y_val;                                                    \
    char typenum;                                                              \
    long n, dim;                                                               \
    void *data, *out;                                                          \
                                                                               \
    /*  Parse the data from Python and try to convert it to a usable format. */\
    if (!PyArg_ParseTuple(args, "Od", &x, &W))                                 \
    {                                                                          \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tspecial_functions.%s\n\n"                            \
                     "\rCould not parse inputs. Legal inputs are:\n"           \
                     "\r\tx:     Numpy Array of real numbers (Floats)\n"       \
                     "\r\tW:    Positive real number (Float)\n\rNotes:\n"      \
                     "\r\tx must be a non-empty one dimensional numpy array.", \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check that W is positive and return with error otherwise.            */\
    if (W <= 0)                                                                \
    {                                                                          \
        PyErr_Format(PyExc_ValueError, "\n\rError Encountered: rss_ringoccs\n" \
                                       "\r\tspecial_functions.%s\n\n"          \
                                       "\rW must be a positive number.\n",     \
                                       VarToString(FuncName));                 \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    if (PyLong_Check(x) || PyFloat_Check(x))                                   \
    {                                                                          \
        x_val = PyFloat_AsDouble(x);                                           \
        y_val = rssringoccs_Double_##CName(x_val, W);                          \
        return PyFloat_FromDouble(y_val);                                      \
    }                                                                          \
    else if (PyList_Check(x))                                                  \
    {                                                                          \
        dim    = PyList_Size(x);                                               \
        output = PyList_New(dim);                                              \
                                                                               \
        for (n=0; n<dim; ++n)                                                  \
        {                                                                      \
            nth_item = PyList_GET_ITEM(x, n);                                  \
            if (!PyFloat_Check(nth_item) && !PyLong_Check(nth_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tspecial_functions.%s\n\n"                    \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            x_val = PyFloat_AsDouble(nth_item);                                \
            y_val = rssringoccs_Double_##CName(x_val, W);                      \
            PyList_SET_ITEM(output, n, PyFloat_FromDouble(y_val));             \
        }                                                                      \
        return output;                                                         \
    }                                                                          \
    else if (!(PyArray_Check(x)))                                              \
    {                                                                          \
        PyErr_Format(PyExc_TypeError,                                          \
             "\n\rError Encountered: rss_ringoccs\n"                           \
             "\r\tspecial_functions.%s\n\n"                                    \
             "\rCould not parse inputs. Legal inputs are:\n"                   \
             "\r\tx: Numpy Array of real numbers (Floats)\n\rNotes:\n"         \
             "\r\tx must be a non-empty one dimensional numpy array.\n",       \
             VarToString(FuncName));                                           \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If you get here, then the input is a numpy array. Grab some useful   */\
    /*  information about the data using numpy's API functions.              */\
    typenum = PyArray_TYPE((PyArrayObject *)x);                                \
    dim     = PyArray_DIMS((PyArrayObject *)x)[0];                             \
    data    = PyArray_DATA((PyArrayObject *)x);                                \
                                                                               \
    /*  Check the inputs to make sure they're valid.                         */\
    if (PyArray_NDIM((PyArrayObject *)x) != 1){                                \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tspecial_functions.%s\n"                              \
                     "\n\rInput is not 1-dimensional.\n",                      \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
    else if (dim == 0){                                                        \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tspecial_functions.%s"                                \
                     "\n\n\rInput numpy array is empty.\n",                    \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    if (typenum == NPY_FLOAT)                                                  \
        out = _get_array_from_two_float(data, W, dim,                          \
                                        rssringoccs_Float_##CName);            \
    else if (typenum == NPY_DOUBLE)                                            \
        out = _get_array_from_two_double(data, W, dim,                         \
                                         rssringoccs_Double_##CName);          \
    else if (typenum == NPY_LONGDOUBLE)                                        \
        out = _get_array_from_two_longdouble(data, W, dim,                     \
                                              rssringoccs_LDouble_##CName);    \
    else                                                                       \
    {                                                                          \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *new_x = PyArray_FromObject(x, NPY_DOUBLE, 1, 1);             \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(new_x))                                                          \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tspecial_functions.%s\n\n"                         \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  If it passed, get a pointer to the data inside the numpy array.  */\
        else                                                                   \
            data = PyArray_DATA((PyArrayObject *)new_x);                       \
                                                                               \
        /*  loop over the data and compute with CName_Double.                */\
        out = _get_array_from_two_double(data, W, dim,                         \
                                         rssringoccs_Double_##CName);          \
        typenum = NPY_DOUBLE;                                                  \
    }                                                                          \
                                                                               \
    output  = PyArray_SimpleNewFromData(1, &dim, typenum, out);                \
    capsule = PyCapsule_New(out, NULL, capsule_cleanup);                       \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.  */   \
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

#define WindowFunctionAlForNumpy(FuncName, CName)                              \
static PyObject * FuncName(PyObject *self, PyObject *args)                     \
{                                                                              \
    PyObject *output, *capsule, *x, *nth_item;                                 \
    double W, alpha, x_val, y_val;                                             \
    char typenum;                                                              \
    long n, dim;                                                               \
    void *data, *out;                                                          \
                                                                               \
    /*  Parse the data from Python and try to convert it to a usable format. */\
    if (!PyArg_ParseTuple(args, "Odd", &x, &W, &alpha))                        \
    {                                                                          \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tdiffrec.special_functions.%s\n\n"                    \
                     "\rCould not parse inputs. Legal inputs are:\n"           \
                     "\r\tx:     Numpy Array of real numbers (Floats)\n"       \
                     "\r\tW:     Positive real number (Float)\n"               \
                     "\r\talpha: Positive real number (Float)"                 \
                     "\rNotes:\n"                                              \
                     "\r\tx must be a non-empty one dimensional numpy array.", \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check that W and alpha are positive.                                 */\
    if (W <= 0)                                                                \
    {                                                                          \
        PyErr_Format(PyExc_ValueError, "\n\rError Encountered: rss_ringoccs\n" \
                                       "\r\tspecial_functions.%s\n\n"          \
                                       "\rW must be a positive number.\n",     \
                                       VarToString(FuncName));                 \
        return NULL;                                                           \
    }                                                                          \
    else if (alpha <= 0)                                                       \
    {                                                                          \
        PyErr_Format(PyExc_ValueError, "\n\rError Encountered: rss_ringoccs\n" \
                                       "\r\tspecial_functions.%s\n\n"          \
                                       "\ralpha must be a positive number.\n", \
                                       VarToString(FuncName));                 \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    if (PyLong_Check(x) || PyFloat_Check(x))                                   \
    {                                                                          \
        x_val = PyFloat_AsDouble(x);                                           \
        y_val = rssringoccs_Double_##CName(x_val, W, alpha);                   \
        return PyFloat_FromDouble(y_val);                                      \
    }                                                                          \
    else if (PyList_Check(x))                                                  \
    {                                                                          \
        dim    = PyList_Size(x);                                               \
        output = PyList_New(dim);                                              \
                                                                               \
        for (n=0; n<dim; ++n)                                                  \
        {                                                                      \
            nth_item = PyList_GET_ITEM(x, n);                                  \
            if (!PyFloat_Check(nth_item) && !PyLong_Check(nth_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tspecial_functions.%s\n\n"                    \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            x_val = PyFloat_AsDouble(nth_item);                                \
            y_val = rssringoccs_Double_##CName(x_val, W, alpha);               \
            PyList_SET_ITEM(output, n, PyFloat_FromDouble(y_val));             \
        }                                                                      \
        return output;                                                         \
    }                                                                          \
    else if (!(PyArray_Check(x)))                                              \
    {                                                                          \
        PyErr_Format(PyExc_TypeError,                                          \
             "\n\rError Encountered: rss_ringoccs\n"                           \
             "\r\tspecial_functions.%s\n\n"                                    \
             "\rCould not parse inputs. Legal inputs are:\n"                   \
             "\r\tx: Numpy Array of real numbers (Floats)\n\rNotes:\n"         \
             "\r\tx must be a non-empty one dimensional numpy array.\n",       \
             VarToString(FuncName));                                           \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If you get here, then the input is a numpy array. Grab some useful   */\
    /*  information about the data using numpy's API functions.              */\
    typenum = PyArray_TYPE((PyArrayObject *)x);                                \
    dim     = PyArray_DIMS((PyArrayObject *)x)[0];                             \
    data    = PyArray_DATA((PyArrayObject *)x);                                \
                                                                               \
    /*  Check the inputs to make sure they're valid.                         */\
    if (PyArray_NDIM((PyArrayObject *)x) != 1){                                \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tspecial_functions.%s\n"                              \
                     "\n\rInput is not 1-dimensional.\n",                      \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
    else if (dim == 0){                                                        \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\rError Encountered: rss_ringoccs\n"                   \
                     "\r\tspecial_functions.%s"                                \
                     "\n\n\rInput numpy array is empty.\n",                    \
                     VarToString(FuncName));                                   \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    if (typenum == NPY_FLOAT)                                                  \
        out = _get_array_from_three_float(data, W, alpha, dim,                 \
                                         rssringoccs_Float_##CName);           \
    else if (typenum == NPY_DOUBLE)                                            \
        out = _get_array_from_three_double(data, W, alpha, dim,                \
                                           rssringoccs_Double_##CName);        \
    else if (typenum == NPY_LONGDOUBLE)                                        \
        out = _get_array_from_three_longdouble(data, W, alpha, dim,            \
                                               rssringoccs_LDouble_##CName);   \
    else                                                                       \
    {                                                                          \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *new_x = PyArray_FromObject(x, NPY_DOUBLE, 1, 1);             \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(new_x))                                                          \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tspecial_functions.%s\n\n"                         \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  If it passed, get a pointer to the data inside the numpy array.  */\
        else                                                                   \
            data = PyArray_DATA((PyArrayObject *)new_x);                       \
                                                                               \
        /*  loop over the data and compute with CName_Double.                */\
        out = _get_array_from_three_double(data, W, alpha,                     \
                                           dim, rssringoccs_Double_##CName);   \
        typenum = NPY_DOUBLE;                                                  \
    }                                                                          \
                                                                               \
    output  = PyArray_SimpleNewFromData(1, &dim, typenum, out);                \
    capsule = PyCapsule_New(out, NULL, capsule_cleanup);                       \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.  */   \
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

#define MinMaxFunctionForNumpy(FuncName, CName)                                \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    PyArrayObject *arr;                                                        \
    PyObject *tuple = PyTuple_GetItem(args, 0);                                \
                                                                               \
    if (PyLong_Check(tuple))                                                   \
    {                                                                          \
        long max;                                                              \
        if (PyArg_ParseTuple(args, "l", &max))return PyLong_FromLong(max);     \
        else                                                                   \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                         "\n\r\trss_ringoccs.diffrec.special_functions.%s\n"   \
                         "\r\t\tCould not parse int type input.",              \
                         VarToString(FuncName));                               \
            return NULL;                                                       \
        }                                                                      \
    }                                                                          \
    else if (PyFloat_Check(tuple))                                             \
    {                                                                          \
        double max;                                                            \
        if(PyArg_ParseTuple(args, "d", &max)) return PyFloat_FromDouble(max);  \
        else                                                                   \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                         "\n\r\trss_ringoccs.diffrec.special_functions.%s\n"   \
                         "\r\t\tCould not parse float type input.",            \
                         VarToString(FuncName));                               \
            return NULL;                                                       \
        }                                                                      \
    }                                                                          \
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr))                \
    {                                                                          \
        long typenum, dim;                                                     \
        void *data;                                                            \
                                                                               \
        /*  Check to make sure input isn't zero dimensional!                 */\
        if (PyArray_NDIM(arr) != 1) goto FAIL;                                 \
                                                                               \
        /*  Useful information about the data.                               */\
        typenum = PyArray_TYPE(arr);                                           \
        dim     = PyArray_DIMS(arr)[0];                                        \
        data    = PyArray_DATA(arr);                                           \
                                                                               \
        if (dim == 0) goto FAIL;                                               \
        else if (typenum == NPY_FLOAT)                                         \
            return PyFloat_FromDouble(CName##_Float((float *)data, dim));      \
        else if (typenum == NPY_DOUBLE)                                        \
            return PyFloat_FromDouble(CName##_Double((double *)data, dim));    \
        else if (typenum == NPY_LONGDOUBLE)                                    \
            return PyFloat_FromDouble(                                         \
                CName##_LDouble((long double *)data, dim)                      \
            );                                                                 \
        else if (typenum == NPY_BYTE)                                          \
            return PyLong_FromLong(CName##_Char((char *)data, dim));           \
        else if (typenum == NPY_UBYTE)                                         \
            return PyLong_FromLong(CName##_UChar((unsigned char *)data, dim)); \
        else if (typenum == NPY_SHORT)                                         \
            return PyLong_FromLong(CName##_Short((short *)data, dim));         \
        else if (typenum == NPY_USHORT)                                        \
            return PyLong_FromLong(CName##_UShort((unsigned short *)data, dim));\
        else if (typenum == NPY_INT)                                           \
            return PyLong_FromLong(CName##_Int((int *)data, dim));             \
        else if (typenum == NPY_UINT)                                          \
            return PyLong_FromLong(CName##_UInt((unsigned int *)data, dim));   \
        else if (typenum == NPY_LONG)                                          \
            return PyLong_FromLong(CName##_Long((long *)data, dim));           \
        else if (typenum == NPY_ULONG)                                         \
            return PyLong_FromLong(CName##_ULong((unsigned long *)data, dim)); \
        else goto FAIL;                                                        \
    }                                                                          \
    else goto FAIL;                                                            \
    FAIL: {                                                                    \
        PyErr_Format(PyExc_TypeError,                                          \
                     "\n\r\trss_ringoccs.diffrec.special_functions.%s\n"       \
                     "\r\t\tInput should be a one dimensional numpy array of\n"\
                     "\r\t\treal numbers, or a float/int number.\n"            \
                     "\r\t\tExample:\n"                                        \
                     "\r\t\t\t>>> import numpy\n"                              \
                     "\r\t\t\t>>> import special_functions\n"                  \
                     "\r\t\t\t>>> x = numpy.random.rand(100)\n"                \
                     "\r\t\t\t>>> y = special_functions.%s(x)\n\n"             \
                     "\r\t\tNOTE:\n"                                           \
                     "\r\t\t\tOnly one dimensional numpy arrays are allowed.\n"\
                     "\r\t\t\tComplex numbers are not allowed. If the input\n" \
                     "\r\t\t\tis a single floating point or integer number,\n" \
                     "\r\t\t\tthe output will simply be that number.",         \
                     VarToString(FuncName), VarToString(FuncName));            \
        return NULL;                                                           \
    };                                                                         \
}


static PyObject *besselI0(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Bessel_I0;
    c_funcs.double_func = rssringoccs_Double_Bessel_I0;
    c_funcs.ldouble_func = rssringoccs_LDouble_Bessel_I0;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = rssringoccs_CDouble_Bessel_I0;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *besselJ0(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Bessel_J0;
    c_funcs.double_func = rssringoccs_Double_Bessel_J0;
    c_funcs.ldouble_func = rssringoccs_LDouble_Bessel_J0;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *sinc(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Sinc;
    c_funcs.double_func = rssringoccs_Double_Sinc;
    c_funcs.ldouble_func = rssringoccs_LDouble_Sinc;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *fresnel_sin(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Fresnel_Sin;
    c_funcs.double_func = rssringoccs_Double_Fresnel_Sin;
    c_funcs.ldouble_func = rssringoccs_LDouble_Fresnel_Sin;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *fresnel_cos(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Fresnel_Cos;
    c_funcs.double_func = rssringoccs_Double_Fresnel_Cos;
    c_funcs.ldouble_func = rssringoccs_LDouble_Fresnel_Cos;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *lambertw(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_LambertW;
    c_funcs.double_func = rssringoccs_Double_LambertW;
    c_funcs.ldouble_func = rssringoccs_LDouble_LambertW;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *wavelength_to_wavenumber(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Wavelength_To_Wavenumber;
    c_funcs.double_func = rssringoccs_Double_Wavelength_To_Wavenumber;
    c_funcs.ldouble_func = rssringoccs_LDouble_Wavelength_To_Wavenumber;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *frequency_to_wavelength(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Frequency_To_Wavelength;
    c_funcs.double_func = rssringoccs_Double_Frequency_To_Wavelength;
    c_funcs.ldouble_func = rssringoccs_LDouble_Frequency_To_Wavelength;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

static PyObject *resolution_inverse(PyObject *self, PyObject *args)
{
    rssringoccs_Generic_Function_Obj c_funcs;

    c_funcs.long_func = NULL;
    c_funcs.float_func = rssringoccs_Float_Resolution_Inverse;
    c_funcs.double_func = rssringoccs_Double_Resolution_Inverse;
    c_funcs.ldouble_func = rssringoccs_LDouble_Resolution_Inverse;
    c_funcs.cdouble_from_real_func = NULL;
    c_funcs.cdouble_from_complex_func = NULL;

    return rssringoccs_Get_Py_Func_From_C(self, args, &c_funcs);
}

WindowFunctionForNumpy(rect,   Rect_Window)
WindowFunctionForNumpy(coss,   Coss_Window)
WindowFunctionForNumpy(kb20,   Kaiser_Bessel_2_0)
WindowFunctionForNumpy(kb25,   Kaiser_Bessel_2_5)
WindowFunctionForNumpy(kb35,   Kaiser_Bessel_3_5)
WindowFunctionForNumpy(kbmd20, Modified_Kaiser_Bessel_2_0)
WindowFunctionForNumpy(kbmd25, Modified_Kaiser_Bessel_2_5)
WindowFunctionForNumpy(kbmd35, Modified_Kaiser_Bessel_3_5)
WindowFunctionAlForNumpy(kbal, Kaiser_Bessel)
WindowFunctionAlForNumpy(kbmdal, Modified_Kaiser_Bessel)
MinMaxFunctionForNumpy(min, rssringoccs_Min)
MinMaxFunctionForNumpy(max, rssringoccs_Max)

static PyObject *compute_norm_eq(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long normeq;
        if (PyArg_ParseTuple(args, "l", &normeq))
            return PyLong_FromLong(normeq);
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double normeq;
        if (PyArg_ParseTuple(args, "d", &normeq)){
            return PyFloat_FromDouble(normeq);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        long typenum, dim;
        void *data;

        /*  Check to make sure input isn't zero dimensional!                  */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "\n\trss_ringoccs.diffrec.special_functions.compute_norm_eq\n"
                "\r\t\tInput must be a one-dimensional array."
            );
            return NULL;
        }

        /* Useful information about the data.                                 */
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);
        if (dim == 0){
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tInput is zero dimensional."
            );
            return NULL;
        }
        if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_Float((float *)data, dim)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_Double((double *)data, dim)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_LDouble((long double *)data, dim)
            );
        }
        else if (typenum == NPY_SHORT){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_Short((short *)data, dim)
            );
        }
        else if (typenum == NPY_INT){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_Int((int *)data, dim)
            );
        }
        else if (typenum == NPY_LONG){
            return PyFloat_FromDouble(
                rssringoccs_Normeq_Long((long *)data, dim)
            );
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tInput should be a numpy array of numbers."
            );
            return NULL;
        }
    }
    else {
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
            "\r\t\tInput should be a numpy array of numbers."
        );
        return NULL;
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      where_greater                                                         *
 *  Purpose:                                                                  *
 *      Given a numpy array arr, and a real number threshold, returns an      *
 *      array of alll indices n such that arr[n] > threshold.                 *
 *  Arguments:                                                                *
 *      arr (Numpy Array):                                                    *
 *          An arbitrary numpy array of real numbers.                         *
 *      threshold (double):                                                   *
 *          A real valued number. This is the value such that, for all n such *
 *          that arr[n] > threshold, the index n will be appended to the      *
 *          returning array.                                                  *
 *  Output:                                                                   *
 *      where_arr:                                                            *
 *          A numpy array of integer corresponding to the indices n such that *
 *          arr[n] > threshold.                                               *
 *  NOTES:                                                                    *
 *      1.) Values that are equal to threshold will not be included.          *
 *      2.) The input array MUST be one-dimensional.                          *
 *      3.) integer and floating point types are allowed, but complex are not.*
 ******************************************************************************/
static PyObject *where_greater(PyObject *self, PyObject *args)
{
    /*  Declare necessary variables.                                          */
    PyObject *output, *capsule;
    PyArrayObject *arr;
    double threshold;
    unsigned long **where;

    /*  The input is a numpy array, proceed. Otherwise spit out an error.     */
    if (PyArg_ParseTuple(args, "O!d", &PyArray_Type, &arr, &threshold)){

        /*  Declare some more necessary variables.                            */
        long typenum, dim;
        void *data;

        /*  Check to make sure input is one dimensional.                      */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "rss_ringoccs.diffrec.special_functions.where_greater\n"
                "\r\tInput must be a one-dimensional array and a real number."
            );
            return NULL;
        }

        /*   Useful information about the data.                               */
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        /*  Compute the index array corresponding to the indices n such that  *
         *  arr[n] > threshold. The returned pointer is created using malloc  *
         *  so we must remember to free it later to avoid memory leaks.       */
        switch(typenum)
        {
            case NPY_BYTE:
                where = rssringoccs_Where_Greater_Char((char *)data, dim, threshold);
                break;
            case NPY_UBYTE:
                where =
                rssringoccs_Where_Greater_UChar((unsigned char *)data, dim, threshold);
                break;
            case NPY_SHORT:
                where = rssringoccs_Where_Greater_Short((short *)data, dim, threshold);
                break;
            case NPY_USHORT:
                where =
                rssringoccs_Where_Greater_UShort((unsigned short *)data, dim, threshold);
                break;
            case NPY_INT:
                where = rssringoccs_Where_Greater_Int((int *)data, dim, threshold);
                break;
            case NPY_UINT:
                where =
                rssringoccs_Where_Greater_UInt((unsigned int *)data, dim, threshold);
                break;
            case NPY_LONG:
                where = rssringoccs_Where_Greater_Long((long *)data, dim, threshold);
                break;
            case NPY_ULONG:
                where =
                rssringoccs_Where_Greater_ULong((unsigned long *)data, dim, threshold);
                break;
            case NPY_FLOAT:
                where = rssringoccs_Where_Greater_Float((float *)data, dim, threshold);
                break;
            case NPY_DOUBLE:
                where = rssringoccs_Where_Greater_Double((double *)data, dim, threshold);
                break;
            case NPY_LONGDOUBLE:
                where =
                rssringoccs_Where_Greater_LDouble((long double *)data, dim, threshold);
                break;
            default:
                PyErr_Format(
                    PyExc_TypeError,
                    "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
                    "\r\t\tInput numpy array should be real valued."
                );
                return NULL;
        }

        /*  The first element of where is the array of indices.               */
        unsigned long *where_arr = where[0];

        /*  The second element is the number of points in the index array.    */
        long where_dim = (long)*where[1];

        /*  Create a Numpy array object to be passed back to Python.          */
        output =
        PyArray_SimpleNewFromData(1, &where_dim, NPY_LONG, (void *)where_arr);

        /*  This frees the variable at the Python level once it's destroyed.  */
        capsule = PyCapsule_New(where_arr, NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)output, capsule);

        /*  Where is created using malloc, so make sure to free it.           */
        free(where);

        /*  Return the results to Python.                                     */
        return Py_BuildValue("N", output);
    }
    else {
        /*  If he input is not a numpy array, return to Python with an error. */
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.special_functions.where_greater\n"
            "\r\t\tInput should be a real numpy array and a real number."
        );
        return NULL;
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      where_lesser                                                          *
 *  Purpose:                                                                  *
 *      Given a numpy array arr, and a real number threshold, returns an      *
 *      array of all indices n such that arr[n] < threshold.                  *
 *  Arguments:                                                                *
 *      arr (Numpy Array):                                                    *
 *          An arbitrary numpy array of real numbers.                         *
 *      threshold (double):                                                   *
 *          A real valued number. This is the value such that, for all n such *
 *          that arr[n] < threshold, the index n will be appended to the      *
 *          returning array.                                                  *
 *  Output:                                                                   *
 *      where_arr:                                                            *
 *          A numpy array of integer corresponding to the indices n such that *
 *          arr[n] < threshold.                                               *
 *  NOTES:                                                                    *
 *      1.) Values that are equal to threshold will not be included.          *
 *      2.) The input array MUST be one-dimensional.                          *
 *      3.) integer and floating point types are allowed, but complex are not.*
 ******************************************************************************/
static PyObject *where_lesser(PyObject *self, PyObject *args)
{
    /*  Declare necessary variables.                                          */
    PyObject *output, *capsule;
    PyArrayObject *arr;
    double threshold;
    unsigned long **where;
    long typenum, dim;

    /*  The input is a numpy array, proceed. Otherwise spit out an error.     */
    if (PyArg_ParseTuple(args, "O!d", &PyArray_Type, &arr, &threshold)){

        /*  Check to make sure input is one dimensional.                      */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "rss_ringoccs.diffrec.special_functions.where_greater\n"
                "\r\tInput must be a one-dimensional array and a real number."
            );
            return NULL;
        }

        /*   Useful information about the data.                               */
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];

        /*  Compute the index array corresponding to the indices n such that  *
         *  arr[n] > threshold. The returned pointer is created using malloc  *
         *  so we must remember to free it later to avoid memory leaks.       */
        if (typenum == NPY_BYTE){
            char *data = (char *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Char(data, dim, threshold);
        }
        else if (typenum == NPY_UBYTE){
            unsigned char *data = (unsigned char *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_UChar(data, dim, threshold);
        }
        else if (typenum == NPY_SHORT){
            short *data = (short *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Short(data, dim, threshold);
        }
        else if (typenum == NPY_USHORT){
            unsigned short *data = (unsigned short *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_UShort(data, dim, threshold);
        }
        else if (typenum == NPY_INT){
            int *data = (int *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Int(data, dim, threshold);
        }
        else if (typenum == NPY_UINT){
            unsigned int *data = (unsigned int *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_UInt(data, dim, threshold);
        }
        else if (typenum == NPY_LONG){
            long *data = (long *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Long(data, dim, threshold);
        }
        else if (typenum == NPY_ULONG){
            unsigned long *data = (unsigned long *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_ULong(data, dim, threshold);
        }
        else if (typenum == NPY_FLOAT){
            float *data = (float *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Float(data, dim, threshold);
        }
        else if (typenum == NPY_DOUBLE){
            double *data = (double *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_Double(data, dim, threshold);
        }
        else if (typenum == NPY_LONGDOUBLE){
            long double *data = (long double *)PyArray_DATA(arr);
            where = rssringoccs_Where_Lesser_LDouble(data, dim, threshold);
        }
        else {
            /*  If he input is not a numpy array, return to Python with an error. */
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
                "\r\t\tInput numpy array should be real valued."
            );
            return NULL;
        }

        /*  The first element of where is the array of indices.               */
        unsigned long *where_arr = where[0];

        /*  The second element is the number of points in the index array.    */
        long where_dim = (long)*where[1];

        /*  Create a Numpy array object to be passed back to Python.          */
        output  = PyArray_SimpleNewFromData(1, &where_dim, NPY_LONG,
                                            (void *)where_arr);

        /*  This frees the variable at the Python level once it's destroyed.  */
        capsule = PyCapsule_New(where_arr, NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)output, capsule);

        /*  Where is created using malloc, so make sure to free it.           */
        free(where);

        /*  Return the results to Python.                                     */
        return Py_BuildValue("N", output);
    }
    else {
        /*  If he input is not a numpy array, return to Python with an error. */
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
            "\r\t\tInput should be a numpy array of numbers and a real number."
        );
        return NULL;
    }
}

static PyObject *window_norm(PyObject *self, PyObject *args){
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);
    double dx, f_scale;

    if (PyLong_Check(tuple) || PyFloat_Check(tuple))
    {
        double ker;
        PyArg_ParseTuple(args, "ddd", &ker, &dx, &f_scale);
        return PyFloat_FromDouble(
            rssringoccs_Double_Window_Normalization(&ker, 1, dx, f_scale)
        );
    }
    else if (PyArg_ParseTuple(args, "O!dd", &PyArray_Type,
                              &arr, &dx, &f_scale)){
        long typenum, dim;
        void *data;

        /* Check to make sure input isn't zero dimensional!                   */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError,
                         "rss_ringoccs.diffrec.math_functions.min\n"
                         "\rInput must be a one-dimensional array.");
            return NULL;
        }

        /* Useful information about the data.                                 */
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (typenum == NPY_CDOUBLE){
            return PyFloat_FromDouble(
                rssringoccs_Complex_Window_Normalization(
                    (rssringoccs_ComplexDouble *)data, dim, dx, f_scale
                )
            );
        }
        else if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                rssringoccs_Float_Window_Normalization((float *)data, dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                rssringoccs_Double_Window_Normalization((double *)data, dim, dx, f_scale)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                rssringoccs_LDouble_Window_Normalization((long double *)data,
                                                dim, dx, f_scale)
            );
        }
        else if ((typenum == NPY_INT)       || (typenum == NPY_UINT)        ||
                 (typenum == NPY_SHORT)     || (typenum == NPY_USHORT)      ||
                 (typenum == NPY_LONG)      || (typenum == NPY_ULONG))
        {
            return PyFloat_FromDouble(
                rssringoccs_Double_Window_Normalization((double *)data, dim, dx, f_scale)
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

static PyMethodDef special_functions_methods[] =
{
    {
        "coss",
        coss,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.coss\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the squared cosine window function.\n\r\t"
        "Arguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for coss(x)."
        "\n\r\t\tW (float):\n\r\t\t\t"
        "The width of the window function.\n\r\t"
        "Outputs:\n\r\t\t"
        "coss (numpy.ndarray):\n\r\t\t\t"
        "The squared cosine function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0\n\r\t\t"
        ">>> y = special_functions.coss(x, W)"
    },
    {
        "rect",
        rect,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.rect\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the rectangular window function.\n\r\t"
        "Arguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for rect(x)\n\r\t"
        "W (float):\n\r\t\t\t"
        "The width of the window function.\n\r\t"
        "Outputs:\n\r\t\t"
        "rect (numpy.ndarray):\n\r\t\t\t"
        "The rect function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.rect(x, W)"
    },
    {
        "kb20",
        kb20,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kb20\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Kaiser-Bessel function with alpha = 2.0."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kb20(x)\n\r\t"
        "W (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kb20 (numpy.ndarray):\n\r\t\t\t"
        "The kb20 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.kb20(x, W)"
    },
    {
        "kb25",
        kb25,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kb25\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Kaiser-Bessel function with alpha = 2.5."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kb25(x)\n\r\t"
        "W (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kb25 (numpy.ndarray):\n\r\t\t\t"
        "The kb25 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.kb25(x, W)"
    },
    {
        "kb35",
        kb35,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kb35\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Kaiser-Bessel function with alpha = 3.5."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kb35(x)\n\r\t"
        "W (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kb35 (numpy.ndarray):\n\r\t\t\t"
        "The kb35 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.kb35(x, W)"
    },
    {
        "kbal",
        kbal,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kbal\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Kaiser-Bessel function with arbitrary alpha."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kbal(x)\n\r\t"
        "W (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kbal (numpy.ndarray):\n\r\t\t\t"
        "The kbal function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0\n\r\t\t"
        ">>> alpha = 1.8\n\r\t\t"
        ">>> y = special_functions.kbal(x, W, alpha)"
    },
    {
        "kbmd20",
        kbmd20,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kbmd20\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Modified Kaiser-Bessel function with alpha = 2.0."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kbmd20(x)"
        "\n\r\tW (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kbmd20 (numpy.ndarray):\n\r\t\t\t"
        "The kbmd20 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0\n\r\t\t"
        ">>> y = special_functions.kbmd20(x, W)"
    },
    {
        "kbmd25",
        kbmd25,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kbmd25\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Modified Kaiser-Bessel function with alpha = 2.5."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kbmd25(x)"
        "\n\r\tW (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kbmd25 (numpy.ndarray):\n\r\t\t\t"
        "The kbmd25 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0\n\r\t\t"
        ">>> y = special_functions.kbmd25(x, W)"
    },
    {
        "kbmd35",
        kbmd35,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kbmd35\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Modified Kaiser-Bessel function with alpha = 3.5."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kbmd35(x)"
        "\n\r\tW (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kbmd35 (numpy.ndarray):\n\r\t\t\t"
        "The kbmd35 function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.kbmd35(x, W)"
    },
    {
        "kbmdal",
        kbmdal,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.kbmdal\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Modified Kaiser-Bessel function with arbitrary alpha."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for kbmdal(x)"
        "\n\r\tW (float):\n\r\t\t\t"
        "The width of the Kaiser-Bessel window.\n\r\t"
        "Outputs:\n\r\t\t"
        "kbmdal (numpy.ndarray):\n\r\t\t\t"
        "The kbmdal function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-20,20,0.1)\n\r\t\t"
        ">>> W = 10.0"
        ">>> y = special_functions.kbmdal(x, W)"
    },
    {
        "besselJ0",
        besselJ0,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.besselJ0\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Zeroth Bessel Function of the First Kind (J0).\n\r\t"
        "Arguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for J_0(x)\n\r\t"
        "Outputs:\n\r\t\t"
        "J0\n\r\t\t\t"
        "The Bessel Function J0 as a function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> y = special_functions.besselJ0(x)"
    },
    {
        "besselI0",
        besselI0,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.besselI0\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Zeroth Modified Bessel Function of the First Kind (I0)."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for I_0(x)\n\r\t"
        "Outputs:\n\r\t\t"
        "I0 (numpy.ndarray):\n\r\t\t\t"
        "The Bessel Function I0 as a function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> y = special_functions.besselI0(x)"
    },
    {
        "fresnel_sin",
        fresnel_sin,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.fresnel_sin\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Fresnel sine function."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for fresnel_sin(x)"
        "\n\r\tOutputs:\n\r\t\t"
        "fresnel_sin (numpy.ndarray):\n\r\t\t\t"
        "The fresnel sine function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0, 10, 0.01)\n\r\t\t"
        ">>> y = special_functions.fresnel_sin(x)"
    },
    {
        "fresnel_cos",
        fresnel_cos,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.fresnel_cos\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Fresnel cosine function."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for fresnel_cos(x)"
        "\n\r\tOutputs:\n\r\t\t"
        "fresnel_cos (numpy.ndarray):\n\r\t\t\t"
        "The fresnel sine function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0, 10, 0.01)\n\r\t\t"
        ">>> y = special_functions.fresnel_cos(x)"
    },
    {
        "lambertw",
        lambertw,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.lambertw\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the Lambert W function, inverse of x*exp(x)."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable.\n\r\t"
        "Outputs:\n\r\t\t"
        "y (numpy.ndarray):\n\r\t\t\t"
        "The Lambert W function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> y = special_functions.lambertw(x)"
    },
    {
        "sinc",
        sinc,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.sinc\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute the sinc function, sin(x)/x."
        "\n\r\tArguments\n\r\t\t"
        "x (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers. Independent variable for sinc(x)\n\r\t"
        "Outputs:\n\r\t\t"
        "sinc (numpy.ndarray):\n\r\t\t\t"
        "The sinc function of x.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> y = special_functions.sinc(x)"
    },
    {
        "compute_norm_eq",
        compute_norm_eq,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.compute_norm_eq\n\r\t"
        "Purpose:\n\r\t\t"
        "Compute normalized equivalenth width of a given function.\n\r\t"
        "Arguments:\n\r\t\t"
        "w_func (*numpy.ndarray*):\n\r\t\t\t"
        "Function to compute the normalized equivalent width.\n\r\t"
        "Outputs:\n\r\t\t"
        "normeq (*float*):\n\r\t\t\t"
        "The normalized equivalent width of w_func.\n\r\t"
        "Notes:\n\r\t\t"
        "The normalized equivalent width is computed using Riemann\n\r\t\t"
        "sums to approximate integrals. Therefore large dx values\n\r\t\t"
        "(Spacing between points) will result in an inaccurate\n\r\t\t"
        "normeq. One should keep this in mind during calculations.\n\r\t"
        "Examples:\n\r\t\t"
        "Compute the Kaiser-Bessel 2.5 window of width 20km and\n\r\t\t"
        "spacing 0.1 and compute the normalized equivalent width:\n\r\t\t\t"
        ">>> import special_functions\n\r\t\t\t"
        ">>> import numpy\n\r\t\t\t"
        ">>> x = numpy.arange(-10, 10, 0.1)\n\r\t\t\t"
        ">>> w = special_functions.kb25(x, 20)\n\r\t\t\t"
        ">>> special_functions.compute_norm_eq(w)\n\r\t\t\t"
        "1.651925635118099\n\r\t\t"
        "In contrast, the actual value is 1.6519208. Compute the\n\r\t\t"
        "normalized equivalent width for the squared cosine window of\n\r\t\t"
        "width 20 and spacing 0.25.\n\r\t\t\t"
        ">>> import special_functions\n\r\t\t\t"
        ">>> import numpy\n\r\t\t\t"
        ">>> x = numpy.arange(-10, 10, 0.25)\n\r\t\t\t"
        ">>> w = special_functions.kb25(x, 20)\n\r\t\t\t"
        ">>> special_functions.compute_norm_eq(w)\n\r\t\t\t"
        "1.5000000000000013\n\r\t\t"
        "The normalized equivalent width of the squared cosine\n\r\t\t"
        "function can be computed exactly using standard methods\n\r\t\t"
        "from a calculus course. It's value is exactly 1.5."
    },
    {
        "frequency_to_wavelength",
        frequency_to_wavelength,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.frequency_to_wavelength\n\r\t"
        "Purpose\n\r\t\t"
        "Convert a non-zero frequency to the equivalent wavelength.\n\r\t"
        "Arguments:\n\r\t\t"
        "frequency (numpy.ndarray/int/float/list):\n\r\t\t\t"
        "A numpy array or list of real numbers, input frequency in Hz.\n\r\t\t"
        "Outputs:\n\r\t\t"
        "wavelength (numpy.ndarry/list/float):\n\r\t\t\t"
        "The corresponding wavelength in km.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(1,10,0.1)"
        ">>> y = special_functions.frequency_to_wavelength(x)"
    },
    {
        "max",
        max,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.max\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the maximum of a numpy array.\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t"
        "Outputs:\n\r\t\t"
        "max (int or float):\n\r\t\t\t"
        "The maximum value of the input array.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.random.rand(100)\n\r\t\t"
        ">>> y = special_functions.max(x)"
    },
    {
        "min",
        min,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.min\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the minimum of a numpy array.\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t"
        "Outputs:\n\r\t\t"
        "min (int or float):\n\r\t\t\t"
        "The minimum value of the input array.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.random.rand(100)\n\r\t\t"
        ">>> y = special_functions.min(x)"
    },
    {
        "wavelength_to_wavenumber",
        wavelength_to_wavenumber,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.wavelength_to_wavenumber\n\r\t"
        "Purpose\n\r\t\t"
        "Convert a non-zero wavelength to the equivalent wavenumber.\n\r\t"
        "Arguments:\n\r\t\t"
        "wavelength (numpy.ndarray/int/float/list):\n\r\t\t\t"
        "A numpy array or list of real numbers, the input wavelength.\n\r\t\t"
        "Outputs:\n\r\t\t"
        "wavenumber (float):\n\r\t\t\t"
        "The corresponding wavenumber.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(1,10,0.1)"
        ">>> y = special_functions.wavelength_to_wavenumber(x)"
    },
    {
        "resolution_inverse",
        resolution_inverse,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.resolution_inverse\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the inverse of y = x/(exp(-x)+x-1).\n\r\t"
        "Arguments:\n\r\t\t"
        "x (numpy.ndarray/int/float/list):\n\r\t\t\t"
        "A numpy array or list of real numbers.\n\r\t\t"
        "Outputs:\n\r\t\t"
        "y (float):\n\r\t\t\t"
        "The inverse of x/(exp(-x)+x-1).\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(1,10,0.1)"
        ">>> y = special_functions.resolution_inverse(x)"
    },
    {
        "where_greater",
        where_greater,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.where_greater\n\r\t"
        "Purpose\n\r\t\t"
        "Given a real-valued numpy array arr, and a real number\n\r\t\t"
        "threshold, compute the indices n such that arr[n] > threshold\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t\t"
        "threshold (int or float):\n\r\t\t\t"
        "The threshold value for comparing arr with."
        "Outputs:\n\r\t\t"
        "where_arr (numpy.ndarray):\n\r\t\t\t"
        "The array of indices such that arr[n] > threshold.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-5, 5, 0.01)\n\r\t\t"
        ">>> y = special_functions.where_greater(x, 1.0)"
    },
    {
        "where_lesser",
        where_lesser,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.where_lesser\n\r\t"
        "Purpose\n\r\t\t"
        "Given a real-valued numpy array arr, and a real number\n\r\t\t"
        "threshold, compute the indices n such that arr[n] < threshold\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t\t"
        "threshold (int or float):\n\r\t\t\t"
        "The threshold value for comparing arr with."
        "Outputs:\n\r\t\t"
        "where_arr (numpy.ndarray):\n\r\t\t\t"
        "The array of indices such that arr[n] < threshold.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(-5, 5, 0.01)\n\r\t\t"
        ">>> y = special_functions.where_lesser(x, 1.0)"
    },
    {
        "window_norm",
        window_norm,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "special_functions.window_norm\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the window normalization scheme.\n\r\t"
        "Arguments:\n\r\t\t"
        "ker (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers, the input function.\n\r\t\t"
        "dx (int or float):\n\r\t\t\t"
        "The sample spacing of the input function.\n\r\t\t"
        "f_scale (int or float):\n\r\t\t\t"
        "The Fresnel scale in the same units as dx.\n\r\t\t"
        "Outputs:\n\r\t\t"
        "window_norm (float):\n\r\t\t\t"
        "The normalization factor.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> dx = 0.1"
        ">>> x = numpy.arange(-10,10,dx)\n\r\t\t"
        ">>> ker = special_functions.coss(x, 5)"
        ">>> f_scale = 0.5"
        ">>> y = special_functions.window_norm(ker, dx, f_scale)"
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "special_functions",
    NULL,
    -1,
    special_functions_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_special_functions(void)
{
    PyObject *m = PyModule_Create(&moduledef);
    if (!m) return NULL;

    import_array();

    return m;
}
