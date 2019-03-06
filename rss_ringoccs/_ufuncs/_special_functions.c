#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

/* Define Coefficients for the Fresnel Sine Taylor Expansion. */
#define FRESNEL_SINE_TAYLOR_00 0.3333333333333333
#define FRESNEL_SINE_TAYLOR_01 -0.023809523809523808
#define FRESNEL_SINE_TAYLOR_02 0.0007575757575757576
#define FRESNEL_SINE_TAYLOR_03 -1.3227513227513228e-05
#define FRESNEL_SINE_TAYLOR_04 1.4503852223150468e-07
#define FRESNEL_SINE_TAYLOR_05 -1.0892221037148573e-09
#define FRESNEL_SINE_TAYLOR_06 5.9477940136376354e-12
#define FRESNEL_SINE_TAYLOR_07 -2.466827010264457e-14
#define FRESNEL_SINE_TAYLOR_08 8.032735012415773e-17
#define FRESNEL_SINE_TAYLOR_09 -2.107855191442136e-19
#define FRESNEL_SINE_TAYLOR_10 4.5518467589282e-22
#define FRESNEL_SINE_TAYLOR_11 -8.230149299214221e-25
#define FRESNEL_SINE_TAYLOR_12 1.2641078988989164e-27
#define FRESNEL_SINE_TAYLOR_13 -1.669761793417372e-30

/* Define Coefficients for the Fresnel Sine Asymptotic Expansion. */
#define FRESNEL_SINE_ASYM_00 -0.5
#define FRESNEL_SINE_ASYM_01 -0.25
#define FRESNEL_SINE_ASYM_02 0.375
#define FRESNEL_SINE_ASYM_03 0.9375

/* Define Coefficients for the Fresnel Cosine Taylor Expansion. */
#define FRESNEL_COSINE_TAYLOR_00 1.0
#define FRESNEL_COSINE_TAYLOR_01 -0.1
#define FRESNEL_COSINE_TAYLOR_02 0.004629629629629629
#define FRESNEL_COSINE_TAYLOR_03 -0.00010683760683760684
#define FRESNEL_COSINE_TAYLOR_04 1.4589169000933706e-06
#define FRESNEL_COSINE_TAYLOR_05 -1.3122532963802806e-08
#define FRESNEL_COSINE_TAYLOR_06 8.35070279514724e-11
#define FRESNEL_COSINE_TAYLOR_07 -3.9554295164585257e-13
#define FRESNEL_COSINE_TAYLOR_08 1.4483264643598138e-15
#define FRESNEL_COSINE_TAYLOR_09 -4.221407288807088e-18
#define FRESNEL_COSINE_TAYLOR_10 1.0025164934907719e-20f
#define FRESNEL_COSINE_TAYLOR_11 -1.977064753877905e-23
#define FRESNEL_COSINE_TAYLOR_12 3.289260349175752e-26
#define FRESNEL_COSINE_TAYLOR_13 -4.6784835155184856e-29

/* Define Coefficients for the Fresnel Coine Asymptotic Expansion. */
#define FRESNEL_COSINE_ASYM_00 0.5
#define FRESNEL_COSINE_ASYM_01 -0.25
#define FRESNEL_COSINE_ASYM_02 -0.375
#define FRESNEL_COSINE_ASYM_03 0.9375

/* Define Miscellaneous Constants. */
#define SQRT_PI_BY_8 0.626657068657750125603941
#define SQRT_PI_BY_2 1.2533141373155001
#define SQRT_2_BY_PI 0.7978845608028654

/* Define Coefficients of Taylor Expansion of sin(pi/2 x) and cos(pi/2 x). */
#define SIN_PI_BY_2_X_A01 1.5707963267948966
#define SIN_PI_BY_2_X_A03 -0.6459640975062462
#define SIN_PI_BY_2_X_A05 0.07969262624616703
#define SIN_PI_BY_2_X_A07 -0.004681754135318687

#define COS_PI_BY_2_X_A00 1.0
#define COS_PI_BY_2_X_A02 -1.2337005501361697
#define COS_PI_BY_2_X_A04 0.253669507901048
#define COS_PI_BY_2_X_A06 -0.020863480763352957
#define COS_PI_BY_2_X_A08 0.0009192602748394263


static PyMethodDef _special_functions_methods[] = {{NULL, NULL, 0, NULL}};
/*-----------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.  *
 * The are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.               *
 *----------------------------------------------------------------------------*/
double Fresnel_Sine_Func(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 9.0){
        double arg_sq = arg*arg;
        if (arg < 1.0){
            sx = arg_sq * FRESNEL_SINE_TAYLOR_03 + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
        else if (arg < 4.0){
            sx = arg_sq * FRESNEL_SINE_TAYLOR_07 + FRESNEL_SINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
        else{
            sx = arg_sq * FRESNEL_SINE_TAYLOR_13 + FRESNEL_SINE_TAYLOR_12;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_11;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_10;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_09;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_08;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_07;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
    }
    else if (arg < 1.0e16) {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg;
        arg *= arg;
        sinarg *= arg;
        cosarg *= FRESNEL_SINE_ASYM_02*arg + FRESNEL_SINE_ASYM_00;
        sinarg *= FRESNEL_SINE_ASYM_03*arg + FRESNEL_SINE_ASYM_01;

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

double Fresnel_Cosine_Func(double x)
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
            sx = arg_sq * FRESNEL_COSINE_TAYLOR_03 + FRESNEL_COSINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_00;
            return sx*x;
        }
        else if (arg < 4.0){
            sx = arg_sq * FRESNEL_COSINE_TAYLOR_07 + FRESNEL_COSINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_00;
            return sx*x;
        }
        else{
            sx = arg_sq * FRESNEL_COSINE_TAYLOR_13 + FRESNEL_COSINE_TAYLOR_12;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_11;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_10;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_09;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_08;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_07;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_COSINE_TAYLOR_00;
            return sx*x;
        }
    }
    else {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        sinarg *= arg;
        arg *= arg;
        cosarg *= arg;
        cosarg *= FRESNEL_COSINE_ASYM_03*arg + FRESNEL_COSINE_ASYM_01;
        sinarg *= FRESNEL_COSINE_ASYM_02*arg + FRESNEL_COSINE_ASYM_00;

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

double complex Square_Well_Diffraction_Solution(double x, double a,
                                                double b, double F)
{
    double arg1 = SQRT_2_BY_PI*(x-a)/F;
    double arg2 = SQRT_2_BY_PI*(x-b)/F;
    double real_part = Fresnel_Cosine_Func(arg1) - Fresnel_Cosine_Func(arg2);
    double imag_part = Fresnel_Sine_Func(arg1) - Fresnel_Sine_Func(arg2);
    double complex result = real_part + imag_part*_Complex_I;
    result *= SQRT_PI_BY_2;

    return F*(1.0+1.0*_Complex_I + (0.5 - 0.5*_Complex_I)*result);
}
/*
double Fast_Sine(double x){

}
*/

/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above  *
 * functions to be called in Python. Numpy arrays, as well as floating point  *
 * and integer valued arguments may then be passed into these functions for   *
 * improvement in performance, as opposed to the routines written purely in   *
 * Python. Successful compiling requires the Numpy and Python header files.   *
 *----------------------------------------------------------------------------*/
static void double_fresnelsin(char **args, npy_intp *dimensions,
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
        *((double *)out1) = Fresnel_Sine_Func(*(double *)in1);
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void double_fresnelcos(char **args, npy_intp *dimensions,
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
        *((double *)out1) = Fresnel_Cosine_Func(*(double *)in1);
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void complex_sqwellsol(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *x = args[0];
    char *a = args[1];
    char *b = args[2];
    char *F = args[3];
    char *out = args[4];
    npy_intp in_step = steps[0];
    npy_intp out_step = steps[4];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double complex*)out) = Square_Well_Diffraction_Solution(
            *(double *)x, *(double *)a, *(double *)b, *(double *)F
        );
        /*END main ufunc computation*/

        x += in_step;
        out += out_step;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_sin_funcs[1] = {&double_fresnelsin};
PyUFuncGenericFunction fresnel_cos_funcs[1] = {&double_fresnelcos};
PyUFuncGenericFunction sqwellsol_funcs[1] = {&complex_sqwellsol};

/* Input and return types for double input and out.. */
static char double_double_types[2] = {NPY_DOUBLE, NPY_DOUBLE};
static void *PyuFunc_data[1] = {NULL};

/* Input and return types for square_well_diffraction. */
static char sqwellsol_types[5] = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                  NPY_DOUBLE, NPY_COMPLEX128};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_special_functions",
    NULL,
    -1,
    _special_functions_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__special_functions(void)
{
    PyObject *fresnel_sin;
    PyObject *fresnel_cos;
    PyObject *square_well_diffraction;
    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_sin = PyUFunc_FromFuncAndData(fresnel_sin_funcs, PyuFunc_data,
                                          double_double_types, 1, 1, 1,
                                          PyUFunc_None, "fresnel_sin",
                                          "fresnel_sin_docstring", 0);

    fresnel_cos = PyUFunc_FromFuncAndData(fresnel_cos_funcs, PyuFunc_data,
                                          double_double_types, 1, 1, 1,
                                          PyUFunc_None, "fresnel_cos",
                                          "fresnel_cos_docstring", 0);

    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);

    Py_DECREF(fresnel_sin);
    Py_DECREF(fresnel_cos);
    Py_DECREF(square_well_diffraction);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *m, *fresnel_sin, *fresnel_cos, *d;

    m = Py_InitModule("__funcs", _special_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_sin = PyUFunc_FromFuncAndData(fresnel_sin_funcs, PyuFunc_data,
                                          double_double_types, 1, 1, 1,
                                          PyUFunc_None, "fresnel_sin",
                                          "fresnel_sin_docstring", 0);

    fresnel_cos = PyUFunc_FromFuncAndData(fresnel_cos_funcs, PyuFunc_data,
                                          double_double_types, 1, 1, 1,
                                          PyUFunc_None, "fresnel_cos",
                                          "fresnel_cos_docstring", 0);

    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);

    Py_DECREF(fresnel_sin);
    Py_DECREF(fresnel_cos);
    Py_DECREF(square_well_diffraction);
}
#endif