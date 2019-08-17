/*  To avoid compiler warnings about deprecated numpy stuff.                 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* cosine and sine are defined here. */
#include <math.h>

/*  complex data types, as well as _Complex_I, are defined here.             */
#include <complex.h>

/*  Window functions defined here.                                           */
#include "__window_functions.h"

/*  Various header files required for the C-Python API to work.              */
#include "../../../include/Python.h"
#include "../../../include/ndarraytypes.h"
#include "../../../include/ufuncobject.h"

static PyMethodDef _window_functions_methods[] = {{NULL, NULL, 0, NULL}};
/*---------------------------DEFINE PYTHON FUNCTIONS-------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above *
 * functions to be called in Python. Numpy arrays, as well as floating point *
 * and integer valued arguments may then be passed into these functions for  *
 * improvement in performance, as opposed to the routines written purely in  *
 * Python. Successful compiling requires the Numpy and Python header files.  *
 *---------------------------------------------------------------------------*/
static void double_rect(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Rect_Window_Double(x[i], W);
    }
}

static void double_coss(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Coss_Window_Double(x[i], W);
    }
}

static void double_kb20(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_0_Double(x[i], W);
    }
}

static void double_kb25(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_5_Double(x[i], W);
    }
}

static void double_kb35(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_3_5_Double(x[i], W);
    }
}

static void double_kbmd20(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Modified_Kaiser_Bessel_2_0_Double(x[i], W);
    }
}

static void double_kbmd25(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Modified_Kaiser_Bessel_2_5_Double(x[i], W);
    }
}

static void double_kbmd35(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Modified_Kaiser_Bessel_3_5_Double(x[i], W);
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction rect_funcs[1]   = {&double_rect};
PyUFuncGenericFunction coss_funcs[1]   = {&double_coss};
PyUFuncGenericFunction kb20_funcs[1]   = {&double_kb20};
PyUFuncGenericFunction kb25_funcs[1]   = {&double_kb25};
PyUFuncGenericFunction kb35_funcs[1]   = {&double_kb35};
PyUFuncGenericFunction kbmd20_funcs[1] = {&double_kbmd20};
PyUFuncGenericFunction kbmd25_funcs[1] = {&double_kbmd25};
PyUFuncGenericFunction kbmd35_funcs[1] = {&double_kbmd35};

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
    PyObject *rect, *coss, *kb25, *kb20, *kb35;
    PyObject *kbmd25, *kbmd20, *kbmd35;
    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    rect = PyUFunc_FromFuncAndData(rect_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "rect", "rect_docstring", 0);

    coss = PyUFunc_FromFuncAndData(coss_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "coss", "coss_docstring", 0);

    kb25 = PyUFunc_FromFuncAndData(kb25_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb25", "kb25_docstring", 0);

    kb20 = PyUFunc_FromFuncAndData(kb20_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb20", "kb20_docstring", 0);

    kb35 = PyUFunc_FromFuncAndData(kb35_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb35", "kb35_docstring", 0);

    kbmd20 = PyUFunc_FromFuncAndData(kbmd20_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd20",
                                     "kbmd20_docstring", 0);

    kbmd25 = PyUFunc_FromFuncAndData(kbmd25_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd25",
                                     "kbmd25_docstring", 0);

    kbmd35 = PyUFunc_FromFuncAndData(kbmd35_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd35",
                                     "kbmd35_docstring", 0);

    d = PyModule_GetDict(m);
    PyDict_SetItemString(d, "kbmd20", kbmd20);
    PyDict_SetItemString(d, "kbmd25", kbmd25);
    PyDict_SetItemString(d, "kbmd35", kbmd20);
    PyDict_SetItemString(d, "rect", rect);
    PyDict_SetItemString(d, "coss", coss);
    PyDict_SetItemString(d, "kb20", kb20);
    PyDict_SetItemString(d, "kb25", kb25);
    PyDict_SetItemString(d, "kb35", kb20);
    Py_DECREF(kbmd20);
    Py_DECREF(kbmd25);
    Py_DECREF(kbmd35);
    Py_DECREF(rect);
    Py_DECREF(coss);
    Py_DECREF(kb20);
    Py_DECREF(kb25);
    Py_DECREF(kb35);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *rect, *coss, *kb25, *kb20, *kb35;
    PyObject *kbmd25, *kbmd20, *kbmd35;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _window_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    rect = PyUFunc_FromFuncAndData(rect_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "rect", "rect_docstring", 0);

    coss = PyUFunc_FromFuncAndData(coss_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "coss", "coss_docstring", 0);

    kb25 = PyUFunc_FromFuncAndData(kb25_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb25", "kb25_docstring", 0);

    kb20 = PyUFunc_FromFuncAndData(kb20_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb20", "kb20_docstring", 0);

    kb35 = PyUFunc_FromFuncAndData(kb35_funcs, PyuFunc_data, ddd_types, 1, 2, 1,
                                   PyUFunc_None, "kb35", "kb35_docstring", 0);

    kbmd20 = PyUFunc_FromFuncAndData(kbmd20_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd20",
                                     "kbmd20_docstring", 0);

    kbmd25 = PyUFunc_FromFuncAndData(kbmd25_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd25",
                                     "kbmd25_docstring", 0);

    kbmd35 = PyUFunc_FromFuncAndData(kbmd35_funcs, PyuFunc_data, ddd_types,
                                     1, 2, 1, PyUFunc_None, "kbmd35",
                                     "kbmd35_docstring", 0);

    d = PyModule_GetDict(m);
    PyDict_SetItemString(d, "kbmd20", kbmd20);
    PyDict_SetItemString(d, "kbmd25", kbmd25);
    PyDict_SetItemString(d, "kbmd35", kbmd20);
    PyDict_SetItemString(d, "rect", rect);
    PyDict_SetItemString(d, "coss", coss);
    PyDict_SetItemString(d, "kb20", kb20);
    PyDict_SetItemString(d, "kb25", kb25);
    PyDict_SetItemString(d, "kb35", kb20);
    Py_DECREF(kbmd20);
    Py_DECREF(kbmd25);
    Py_DECREF(kbmd35);
    Py_DECREF(rect);
    Py_DECREF(coss);
    Py_DECREF(kb20);
    Py_DECREF(kb25);
    Py_DECREF(kb35);
}
#endif