/*  To avoid compiler warnings about deprecated numpy stuff.        */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* Include the special functions defined in various header files.   */
#include "__fresnel_integrals.h"

/*  Various header files required for the C-Python API to work.     */
#include "../../../include/Python.h"
#include "../../../include/ndarraytypes.h"
#include "../../../include/ufuncobject.h"

static PyMethodDef _fresnel_integrals_methods[] = {{NULL, NULL, 0, NULL}};
/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above  *
 * functions to be called in Python. Numpy arrays, as well as floating point  *
 * and integer valued arguments may then be passed into these functions for   *
 * improvement in performance, as opposed to the routines written purely in   *
 * Python. Successful compiling requires the Numpy and Python header files.   *
 *----------------------------------------------------------------------------*/

/*------------------------------Fresnel Sine----------------------------------*/
static void fresnelsin_taylor2asymp(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_Taylor_to_Asymptotic_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelsin_while(char **args, npy_intp *dimensions,
                             npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_While_to_Asymptotic_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelsin_healdepsthree(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_Heald_Rational_EPS_Minus_Three_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelsin_healdepsfour(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_Heald_Rational_EPS_Minus_Four_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelsin_healdepssix(char **args, npy_intp *dimensions,
                                   npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_Heald_Rational_EPS_Minus_Six_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelsin_healdepseight(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

/*-----------------------------Fresnel Cosine---------------------------------*/
static void fresnelcos_taylor2asmyp(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Cosine_Taylor_to_Asymptotic_Func is defined  *
         *  in _fresnel_cosine.h. Make sure this is in the current directory! */
        *((double *)out) = Fresnel_Cosine_Taylor_to_Asymptotic_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelcos_while(char **args, npy_intp *dimensions,
                             npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Cosine_While_to_Asymptotic_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelcos_healdepsthree(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelcos_healdepsfour(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelcos_healdepssix(char **args, npy_intp *dimensions,
                                   npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnelcos_healdepseight(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double *)out) = Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

/*-------------------------Complex Fresnel Integral---------------------------*/
static void fresnel_taylor2asmyp(char **args, npy_intp *dimensions,
                                 npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Cosine_Taylor_to_Asymptotic_Func is defined  *
         *  in _fresnel_cosine.h. Make sure this is in the current directory! */
        *((double complex *)out) = Fresnel_Taylor_to_Asymptotic_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnel_healdepsthree(char **args, npy_intp *dimensions,
                                  npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double complex *)out) = Fresnel_Heald_Rational_EPS_Minus_Three_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnel_healdepsfour(char **args, npy_intp *dimensions,
                                 npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double complex *)out) = Fresnel_Heald_Rational_EPS_Minus_Four_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnel_healdepssix(char **args, npy_intp *dimensions,
                                npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double complex *)out) = Fresnel_Heald_Rational_EPS_Minus_Six_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

static void fresnel_healdepseight(char **args, npy_intp *dimensions,
                                  npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n          = dimensions[0];
    char *in            = args[0];
    char *out           = args[1];
    npy_intp in_step    = steps[0];
    npy_intp out_step   = steps[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Sine_Taylor_to_Asymptotic_Func is defined in *
         *  _fresnel_sin.h. Make sure this is in the current directory!       */
        *((double complex *)out) = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(
            *(double *)in
        );

        /* Push the pointers forward by the appropriate increment.            */
        in  += in_step;
        out += out_step;
    }
}

/*--------------------------Square Well Diffraction---------------------------*/
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
        *((double complex*)out) = Square_Well_Diffraction_Solution_Func(
            *(double *)x, *(double *)a, *(double *)b, *(double *)F
        );
        /*END main ufunc computation*/

        x += in_step;
        out += out_step;
    }
}

static void complex_invsqwellsol(char **args, npy_intp *dimensions,
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
        *((double complex*)out) = Inverted_Square_Well_Diffraction_Solution_Func(
            *(double *)x, *(double *)a, *(double *)b, *(double *)F
        );
        /*END main ufunc computation*/

        x += in_step;
        out += out_step;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_sin_t2a_funcs[1] = {&fresnelsin_taylor2asymp};
PyUFuncGenericFunction fresnel_sin_while_funcs[1] = {&fresnelsin_while};
PyUFuncGenericFunction fresnel_sin_healde3_funcs[1] = {&fresnelsin_healdepsthree};
PyUFuncGenericFunction fresnel_sin_healde4_funcs[1] = {&fresnelsin_healdepsfour};
PyUFuncGenericFunction fresnel_sin_healde6_funcs[1] = {&fresnelsin_healdepssix};
PyUFuncGenericFunction fresnel_sin_healde8_funcs[1] = {&fresnelsin_healdepseight};

PyUFuncGenericFunction fresnel_cos_t2a_funcs[1] = {&fresnelcos_taylor2asmyp};
PyUFuncGenericFunction fresnel_cos_while_funcs[1] = {&fresnelcos_while};
PyUFuncGenericFunction fresnel_cos_healde3_funcs[1] = {&fresnelcos_healdepsthree};
PyUFuncGenericFunction fresnel_cos_healde4_funcs[1] = {&fresnelcos_healdepsfour};
PyUFuncGenericFunction fresnel_cos_healde6_funcs[1] = {&fresnelcos_healdepssix};
PyUFuncGenericFunction fresnel_cos_healde8_funcs[1] = {&fresnelcos_healdepseight};

PyUFuncGenericFunction fresnel_t2a_funcs[1] = {&fresnel_taylor2asmyp};
PyUFuncGenericFunction fresnel_healde3_funcs[1] = {&fresnel_healdepsthree};
PyUFuncGenericFunction fresnel_healde4_funcs[1] = {&fresnel_healdepsfour};
PyUFuncGenericFunction fresnel_healde6_funcs[1] = {&fresnel_healdepssix};
PyUFuncGenericFunction fresnel_healde8_funcs[1] = {&fresnel_healdepseight};

PyUFuncGenericFunction sqwellsol_funcs[1]       = {&complex_sqwellsol};
PyUFuncGenericFunction invsqwellsol_funcs[1]    = {&complex_invsqwellsol};

/* Input and return types for double input and out.. */
static char double_double_types[2]  = {NPY_DOUBLE, NPY_DOUBLE};
static char double_complex_types[2] = {NPY_DOUBLE, NPY_COMPLEX128};

/* Input and return types for square_well_diffraction. */
static char sqwellsol_types[5] = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                  NPY_DOUBLE, NPY_COMPLEX128};

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
    PyObject *fresnel_sin_taylor_to_asymptotic;
    PyObject *fresnel_sin_while_to_asymptotic;
    PyObject *fresnel_sin_heald_eps_minus_three;
    PyObject *fresnel_sin_heald_eps_minus_four;
    PyObject *fresnel_sin_heald_eps_minus_six;
    PyObject *fresnel_sin_heald_eps_minus_eight;

    PyObject *fresnel_cos_taylor_to_asymptotic;
    PyObject *fresnel_cos_while_to_asymptotic;
    PyObject *fresnel_cos_heald_eps_minus_three;
    PyObject *fresnel_cos_heald_eps_minus_four;
    PyObject *fresnel_cos_heald_eps_minus_six;
    PyObject *fresnel_cos_heald_eps_minus_eight;

    PyObject *fresnel_integral_taylor_to_asymptotic;
    PyObject *fresnel_integral_heald_eps_minus_three;
    PyObject *fresnel_integral_heald_eps_minus_four;
    PyObject *fresnel_integral_heald_eps_minus_six;
    PyObject *fresnel_integral_heald_eps_minus_eight;

    PyObject *square_well_diffraction;
    PyObject *inverse_square_well_diffraction;
    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    /*  Fresnel Sine Functions.     */
    fresnel_sin_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_t2a_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_taylor_to_asymptotic",
        "fresnel_sin_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_sin_while_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_while_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_while_to_asymptotic",
        "fresnel_sin_while_to_asymptotic_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_sin_healde3_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_three",
        "fresnel_sin_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_sin_healde4_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_four",
        "fresnel_sin_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_sin_healde6_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_six",
        "fresnel_sin_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_sin_healde8_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_eight",
        "fresnel_sin_heald_eps_minus_eight_docstring", 
        0
    );

    /*  Fresnel Cosine Functions.   */
    fresnel_cos_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_cos_t2a_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_taylor_to_asymptotic",
        "fresnel_cos_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_cos_while_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_cos_while_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_while_to_asymptotic",
        "fresnel_cos_while_to_asymptotic_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_cos_healde3_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_three",
        "fresnel_cos_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_cos_healde4_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_four",
        "fresnel_cos_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_cos_healde6_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_six",
        "fresnel_cos_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_cos_healde8_funcs,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_eight",
        "fresnel_cos_heald_eps_minus_eight_docstring", 
        0
    );

    /*  Fresnel Integral Functions.     */
    fresnel_integral_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_t2a_funcs,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_taylor_to_asymptotic",
        "fresnel_integral_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_healde3_funcs,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_heald_eps_minus_three",
        "fresnel_integral_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_healde4_funcs,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_heald_eps_minus_four",
        "fresnel_integral_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_sin_healde6_funcs,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_complex_heald_eps_minus_six",
        "fresnel_complex_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_healde8_funcs,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_complex_heald_eps_minus_eight",
        "fresnel_complex_heald_eps_minus_eight_docstring", 
        0
    );


    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    inverse_square_well_diffraction = PyUFunc_FromFuncAndData(
        invsqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "inverse_square_well_diffraction", 
        "inverse_square_well_diffraction_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_taylor_to_asymptotic",
                             fresnel_sin_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_while_to_asymptotic",
                             fresnel_sin_while_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_three",
                             fresnel_sin_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_four",
                             fresnel_sin_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_six",
                             fresnel_sin_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_eight",
                             fresnel_sin_heald_eps_minus_eight);

    PyDict_SetItemString(d, "fresnel_cos_taylor_to_asymptotic",
                             fresnel_cos_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_cos_while_to_asymptotic",
                             fresnel_cos_while_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_three",
                             fresnel_cos_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_four",
                             fresnel_cos_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_six",
                             fresnel_cos_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_eight",
                             fresnel_cos_heald_eps_minus_eight);

    PyDict_SetItemString(d, "fresnel_integral_taylor_to_asymptotic",
                             fresnel_integral_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_three",
                             fresnel_integral_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_four",
                             fresnel_integral_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_six",
                             fresnel_integral_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_eight",
                             fresnel_integral_heald_eps_minus_eight);

    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);


    Py_DECREF(fresnel_sin_taylor_to_asymptotic);
    Py_DECREF(fresnel_sin_while_to_asymptotic);
    Py_DECREF(fresnel_sin_heald_eps_minus_three);
    Py_DECREF(fresnel_sin_heald_eps_minus_four);
    Py_DECREF(fresnel_sin_heald_eps_minus_six);
    Py_DECREF(fresnel_sin_heald_eps_minus_eight);

    Py_DECREF(fresnel_cos_taylor_to_asymptotic);
    Py_DECREF(fresnel_cos_while_to_asymptotic);
    Py_DECREF(fresnel_cos_heald_eps_minus_three);
    Py_DECREF(fresnel_cos_heald_eps_minus_four);
    Py_DECREF(fresnel_cos_heald_eps_minus_six);
    Py_DECREF(fresnel_cos_heald_eps_minus_eight);

    Py_DECREF(fresnel_integral_taylor_to_asymptotic);
    Py_DECREF(fresnel_integral_heald_eps_minus_three);
    Py_DECREF(fresnel_integral_heald_eps_minus_four);
    Py_DECREF(fresnel_integral_heald_eps_minus_six);
    Py_DECREF(fresnel_integral_heald_eps_minus_eight);

    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_sin_taylor_to_asymptotic;
    PyObject *fresnel_sin_while_to_asymptotic;
    PyObject *fresnel_sin_heald_eps_minus_three;
    PyObject *fresnel_sin_heald_eps_minus_four;
    PyObject *fresnel_sin_heald_eps_minus_six;
    PyObject *fresnel_sin_heald_eps_minus_eight;

    PyObject *fresnel_cos_taylor_to_asymptotic;
    PyObject *fresnel_cos_while_to_asymptotic;
    PyObject *fresnel_cos_heald_eps_minus_three;
    PyObject *fresnel_cos_heald_eps_minus_four;
    PyObject *fresnel_cos_heald_eps_minus_six;
    PyObject *fresnel_cos_heald_eps_minus_eight;

    PyObject *fresnel_integral_taylor_to_asymptotic;
    PyObject *fresnel_integral_heald_eps_minus_three;
    PyObject *fresnel_integral_heald_eps_minus_four;
    PyObject *fresnel_integral_heald_eps_minus_six;
    PyObject *fresnel_integral_heald_eps_minus_eight;

    PyObject *square_well_diffraction;
    PyObject *inverse_square_well_diffraction;

    PyObject *m, *d;

    m = Py_InitModule("__funcs", _fresnel_integrals_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

 /*  Fresnel Sine Functions.     */
    fresnel_sin_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_t2a,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_taylor_to_asymptotic",
        "fresnel_sin_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_sin_while_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_while,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_while_to_asymptotic",
        "fresnel_sin_while_to_asymptotic_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_sin_healde3,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_three",
        "fresnel_sin_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_sin_healde4,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_four",
        "fresnel_sin_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_sin_healde6,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_six",
        "fresnel_sin_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_sin_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_sin_healde8,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_sin_heald_eps_minus_eight",
        "fresnel_sin_heald_eps_minus_eight_docstring", 
        0
    );

    /*  Fresnel Cosine Functions.   */
    fresnel_cos_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_cos_t2a,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_taylor_to_asymptotic",
        "fresnel_cos_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_cos_while_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_cos_while,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_while_to_asymptotic",
        "fresnel_cos_while_to_asymptotic_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_cos_healde3,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_three",
        "fresnel_cos_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_cos_healde4,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_four",
        "fresnel_cos_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_cos_healde6,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_six",
        "fresnel_cos_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_cos_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_cos_healde8,
        PyuFunc_data,
        double_double_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_cos_heald_eps_minus_eight",
        "fresnel_cos_heald_eps_minus_eight_docstring", 
        0
    );

    /*  Fresnel Integral Functions.     */
    fresnel_integral_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_t2a,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_taylor_to_asymptotic",
        "fresnel_integral_taylor_to_asymptotic_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_three = PyUFunc_FromFuncAndData(
        fresnel_healde3,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_heald_eps_minus_three",
        "fresnel_integral_heald_eps_minus_three_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_four = PyUFunc_FromFuncAndData(
        fresnel_healde4,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_integral_heald_eps_minus_four",
        "fresnel_integral_heald_eps_minus_four_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_six = PyUFunc_FromFuncAndData(
        fresnel_sin_healde6,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_complex_heald_eps_minus_six",
        "fresnel_complex_heald_eps_minus_six_docstring", 
        0
    );

    fresnel_integral_heald_eps_minus_eight = PyUFunc_FromFuncAndData(
        fresnel_healde8,
        PyuFunc_data,
        double_complex_types,
        1, 1, 1,
        PyUFunc_None,
        "fresnel_complex_heald_eps_minus_eight",
        "fresnel_complex_heald_eps_minus_eight_docstring", 
        0
    );


    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );
    inverse_square_well_diffraction = PyUFunc_FromFuncAndData(
        invsqwellsol_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "inverse_square_well_diffraction", 
        "inverse_square_well_diffraction_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_taylor_to_asymptotic",
                             fresnel_sin_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_while_to_asymptotic",
                             fresnel_sin_while_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_three",
                             fresnel_sin_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_four",
                             fresnel_sin_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_six",
                             fresnel_sin_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_minus_eight",
                             fresnel_sin_heald_eps_minus_eight);

    PyDict_SetItemString(d, "fresnel_cos_taylor_to_asymptotic",
                             fresnel_cos_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_cos_while_to_asymptotic",
                             fresnel_cos_while_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_three",
                             fresnel_cos_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_four",
                             fresnel_cos_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_six",
                             fresnel_cos_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_cos_heald_eps_minus_eight",
                             fresnel_cos_heald_eps_minus_eight);

    PyDict_SetItemString(d, "fresnel_integral_taylor_to_asymptotic",
                             fresnel_integral_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_three",
                             fresnel_integral_heald_eps_minus_three);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_four",
                             fresnel_integral_heald_eps_minus_four);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_six",
                             fresnel_integral_heald_eps_minus_six);
    PyDict_SetItemString(d, "fresnel_integral_heald_eps_minus_eight",
                             fresnel_integral_heald_eps_minus_eight);

    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);


    Py_DECREF(fresnel_sin_taylor_to_asymptotic);
    Py_DECREF(fresnel_sin_while_to_asymptotic);
    Py_DECREF(fresnel_sin_heald_eps_minus_three);
    Py_DECREF(fresnel_sin_heald_eps_minus_four);
    Py_DECREF(fresnel_sin_heald_eps_minus_six);
    Py_DECREF(fresnel_sin_heald_eps_minus_eight);

    Py_DECREF(fresnel_cos_taylor_to_asymptotic);
    Py_DECREF(fresnel_cos_while_to_asymptotic);
    Py_DECREF(fresnel_cos_heald_eps_minus_three);
    Py_DECREF(fresnel_cos_heald_eps_minus_four);
    Py_DECREF(fresnel_cos_heald_eps_minus_six);
    Py_DECREF(fresnel_cos_heald_eps_minus_eight);

    Py_DECREF(fresnel_integral_taylor_to_asymptotic);
    Py_DECREF(fresnel_integral_heald_eps_minus_three);
    Py_DECREF(fresnel_integral_heald_eps_minus_four);
    Py_DECREF(fresnel_integral_heald_eps_minus_six);
    Py_DECREF(fresnel_integral_heald_eps_minus_eight);

    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
}
#endif