/*  To avoid compiler warnings about deprecated numpy stuff.        */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* cosine and sine are defined here. */
#include <math.h>

/*  complex data types, as well as _Complex_I, are defined here.    */
#include <complex.h>

/* Include the Kaiser-Bessel functions.                             */
#include "__kaiser_bessel.h"

/*  Various header files required for the C-Python API to work.     */
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

static PyMethodDef _window_functions_methods[] = {{NULL, NULL, 0, NULL}};
/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above  *
 * functions to be called in Python. Numpy arrays, as well as floating point  *
 * and integer valued arguments may then be passed into these functions for   *
 * improvement in performance, as opposed to the routines written purely in   *
 * Python. Successful compiling requires the Numpy and Python header files.   *
 *----------------------------------------------------------------------------*/
static void double_kb20(char **args, npy_intp *dimensions,
                        npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *x = args[0];
    char *W = args[1];
    char *out = args[2];

    npy_intp x_step = steps[0];
    npy_intp W_step = steps[1];
    npy_intp out_step = steps[2];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Kaiser_Bessel_Window_2_0(*(double *)x, *(double *)W);

        /* Push the pointers forward by the appropriate increment.            */
        x   += x_step;
        W   += W_step;
        out += out_step;
    }
}

static void double_kb25(char **args, npy_intp *dimensions,
                        npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *x = args[0];
    char *W = args[1];
    char *out = args[2];

    npy_intp x_step = steps[0];
    npy_intp W_step = steps[1];
    npy_intp out_step = steps[2];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Kaiser_Bessel_Window_2_5(*(double *)x, *(double *)W);

        /* Push the pointers forward by the appropriate increment.            */
        x   += x_step;
        W   += W_step;
        out += out_step;
    }
}

static void double_kbmd20(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *x = args[0];
    char *W = args[1];
    char *out = args[2];

    npy_intp x_step = steps[0];
    npy_intp W_step = steps[1];
    npy_intp out_step = steps[2];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Modified_Kaiser_Bessel_Window_2_0(
            *(double *)x, *(double *)W
        );

        /* Push the pointers forward by the appropriate increment.            */
        x   += x_step;
        W   += W_step;
        out += out_step;
    }
}

static void double_kbmd25(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *x = args[0];
    char *W = args[1];
    char *out = args[2];

    npy_intp x_step = steps[0];
    npy_intp W_step = steps[1];
    npy_intp out_step = steps[2];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Modified_Kaiser_Bessel_Window_2_5(
            *(double *)x, *(double *)W
        );

        /* Push the pointers forward by the appropriate increment.            */
        x   += x_step;
        W   += W_step;
        out += out_step;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction kb25_funcs[1] = {&double_kb25};
PyUFuncGenericFunction kb20_funcs[1] = {&double_kb20};
PyUFuncGenericFunction kbmd25_funcs[1] = {&double_kbmd25};
PyUFuncGenericFunction kbmd20_funcs[1] = {&double_kbmd20};

/* Input and return types for double input and out.. */
static char ddd_types[3]     = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static void *PyuFunc_data[1] = {NULL};


#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_window_functions",
    NULL,
    -1,
    _window_functions_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__window_functions(void)
{
    PyObject *kb25, *kb20;
    PyObject *kbmd25, *kbmd20;
    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    kb25 = PyUFunc_FromFuncAndData(kb25_funcs, PyuFunc_data, ddd_types,
                                   1, 2, 1, PyUFunc_None, "kb25",
                                   "kb25_docstring", 0);

    kb20 = PyUFunc_FromFuncAndData(kb20_funcs, PyuFunc_data, ddd_types,
                                   1, 2, 1, PyUFunc_None, "kb20",
                                   "kb20_docstring", 0);

    kbmd25 = PyUFunc_FromFuncAndData(kbmd25_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd25",
                                     "kbmd25_docstring", 0);
    kbmd20 = PyUFunc_FromFuncAndData(kbmd20_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd20",
                                     "kbmd20_docstring", 0);


    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "kbmd25", kbmd25);
    PyDict_SetItemString(d, "kbmd20", kbmd20);
    PyDict_SetItemString(d, "kb25", kb25);
    PyDict_SetItemString(d, "kb20", kb20);
    Py_DECREF(kbmd20);
    Py_DECREF(kbmd25);
    Py_DECREF(kb20);
    Py_DECREF(kb25);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *kbmd25, *kbmd20;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _window_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    kbmd25 = PyUFunc_FromFuncAndData(kbmd25_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd25",
                                     "kbmd25_docstring", 0);
    kbmd20 = PyUFunc_FromFuncAndData(kbmd20_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd20",
                                     "kbmd20_docstring", 0);


    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "kbmd25", kbmd25);
    PyDict_SetItemString(d, "kbmd20", kbmd20);
    Py_DECREF(kbmd20);
    Py_DECREF(kbmd25);
}
#endif