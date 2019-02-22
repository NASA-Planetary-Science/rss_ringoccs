
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

double TAYLOR_COEFFICIENTS[20] = {
    0.3333333333333333, -0.023809523809523808,
    0.0007575757575757576, -1.3227513227513228e-05,
    1.4503852223150468e-07, -1.0892221037148573e-09,
    5.9477940136376354e-12, -2.466827010264457e-14,
    8.032735012415773e-17, -2.107855191442136e-19,
    4.5518467589282e-22, -8.230149299214221e-25,
    1.2641078988989164e-27, -1.669761793417372e-30,
    1.9169428621097826e-33, -1.9303572088151077e-36,
    1.7188560628017835e-39, -1.3630412617791397e-42,
    9.687280238870761e-46, -6.205657919637397e-49
};

static PyMethodDef _fresnel_integrals_methods[] = {{NULL, NULL, 0, NULL}};
double SQRT_PI_BY_8 = 0.626657068657750125603941;

/*-----------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.  *
 * The are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.               *
 *----------------------------------------------------------------------------*/
double Fresnel_Sine_Taylor_to_Asymptotic(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx, arg;
    arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        double arg_sq = arg*arg;
        if (arg < 1.0){
            sx = arg_sq * -1.3227513227513228e-05 + 0.0007575757575757576;
            sx = sx*arg_sq - 0.023809523809523808;
            sx = sx*arg_sq + 0.3333333333333333;
            sx *= arg;
            return sx*x;
        }
        else if (arg < 4.0){
            sx = arg_sq * -2.466827010264457e-14 + 5.9477940136376354e-12;
            sx = arg_sq * sx - 1.0892221037148573e-09;
            sx = arg_sq * sx + 1.4503852223150468e-07;
            sx = arg_sq * sx - 1.3227513227513228e-05;
            sx = arg_sq * sx + 0.0007575757575757576;
            sx = arg_sq * sx - 0.023809523809523808;
            sx = arg_sq * sx + 0.3333333333333333;
            sx *= arg;
            return sx*x;
        }
        else{
            sx = arg_sq * -1.669761793417372e-30 + 1.2641078988989164e-27;
            sx = arg_sq * sx - 8.230149299214221e-25;
            sx = arg_sq * sx + 4.5518467589282e-22;
            sx = arg_sq * sx - 2.107855191442136e-19;
            sx = arg_sq * sx + 8.032735012415773e-17;
            sx = arg_sq * sx - 2.466827010264457e-14;
            sx = arg_sq * sx + 5.9477940136376354e-12;
            sx = arg_sq * sx - 1.0892221037148573e-09;
            sx = arg_sq * sx + 1.4503852223150468e-07;
            sx = arg_sq * sx - 1.3227513227513228e-05;
            sx = arg_sq * sx + 0.0007575757575757576;
            sx = arg_sq * sx - 0.023809523809523808;
            sx = arg_sq * sx + 0.3333333333333333;
            sx *= arg;
            return sx*x;
        }
    }
    else {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg*(0.375*arg*arg - 0.5);
        sinarg *= arg*arg*(0.9375*arg*arg - 0.25);

        sx = cosarg + sinarg;
        sx *= x;
        if (x > 0){
            return sx+SQRT_PI_BY_8;
        }
        else {
            return sx-SQRT_PI_BY_8;
        }
    }
}

double Fresnel_Sine_While_to_Asymptotic(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx, arg, x4, EPS;
    arg = x*x;
    x4 = arg*arg;
    EPS = 1.0e-8;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        int i = 0;
        double term = arg*TAYLOR_COEFFICIENTS[0];
        sx = term;
        while (term > EPS){
            /* Odd terms have a negative coefficients.
               Compute two terms at a time to compare with EPS. */
            i += 1;
            arg *= x4;
            term = arg*TAYLOR_COEFFICIENTS[i];
            sx += term;

            // Compute even term.
            i += 1;
            arg *= x4;
            term = arg*TAYLOR_COEFFICIENTS[i];
            sx += term;
        }
        return x*sx;
    }
    else {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        x4 = 1.0/x4;
        cosarg *= -arg*0.5;
        sinarg *= x4*0.25;

        sx = cosarg + sinarg;
        sx *= x;
        if (x > 0){
            return sx+SQRT_PI_BY_8;
        }
        else {
            return sx-SQRT_PI_BY_8;
        }
    }
}

/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above  *
 * functions to be called in Python. Numpy arrays, as well as floating point  *
 * and integer valued arguments may then be passed into these functions for   *
 * improvement in performance, as opposed to the routines written purely in   *
 * Python. Successful compiling requires the Numpy and Python header files.   *
 *----------------------------------------------------------------------------*/
static void double_fresnelsin_1(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_Taylor_to_Asymptotic(*(double *)in1);
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void double_fresnelsin_2(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_While_to_Asymptotic(*(double *)in1);
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_sin_funcs_1[1] = {&double_fresnelsin_1};
PyUFuncGenericFunction fresnel_sin_funcs_2[1] = {&double_fresnelsin_2};

/* Input and return types for double input and out.. */
static char double_double_types[2] = {NPY_DOUBLE, NPY_DOUBLE};
static void *PyuFunc_data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fresnel_integrals",
    NULL,
    -1,
    _fresnel_integrals_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__fresnel_integrals(void)
{
    PyObject *fresnel_sin_1;
    PyObject *fresnel_sin_2;
    PyObject *fresnel_sin_3;
    PyObject *fresnel_sin_4;
    PyObject *fresnel_sin_5;

    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_sin_1 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_1, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_1",
                                            "fresnel_sin_1_docstring", 0);

    fresnel_sin_2 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_2, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_2",
                                            "fresnel_sin_2_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_1", fresnel_sin_1);
    PyDict_SetItemString(d, "fresnel_sin_2", fresnel_sin_2);
    Py_DECREF(fresnel_sin_1);
    Py_DECREF(fresnel_sin_2);
    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_sin_1;
    PyObject *fresnel_sin_2;
    PyObject *fresnel_sin_3;
    PyObject *fresnel_sin_4;
    PyObject *fresnel_sin_5;

    PyObject *m, *d;

    m = Py_InitModule("__funcs", _fresnel_integrals_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_sin_1 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_1, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_1",
                                            "fresnel_sin_1_docstring", 0);

    fresnel_sin_2 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_2, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_2",
                                            "fresnel_sin_2_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_1", fresnel_sin_1);
    PyDict_SetItemString(d, "fresnel_sin_2", fresnel_sin_2);
    Py_DECREF(fresnel_sin_1);
    Py_DECREF(fresnel_sin_2);
}
#endif