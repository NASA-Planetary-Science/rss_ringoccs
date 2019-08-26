/*  To avoid compiler warnings about deprecated numpy stuff.                 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* cosine and sine are defined here. */
#include <math.h>

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>

/*  Diffraction modeling functions, using Fresnel approximation, found here.  */
#include "__fresnel_diffraction.h"

/*  Various header files required for the C-Python API to work.               */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

static PyMethodDef _fresnel_diffraction_methods[] =
{
    {NULL, NULL, 0, NULL}
};

/*************Square Well Diffraction Using Fresnel Approximation**************/

/******************************************************************************
 *  Function:                                                                 *
 *      complex_float_square_well                                             *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      square well, assuming the Fresnel approximation is valid.             *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) This is a wrapper for Square_Well_Diffraction_Solution_Float,     *
 *          which is defined in __fresnel_diffraction.h. This allows Python   *
 *          to that function, and allows for numpy arrays to be passed in. Ti *
 *          relies on the Numpy UFUNC API, as well as the C-Python API.       *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void complex_float_square_well(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    float *x  =  (float *)args[0];
    float  a  = *(float *)args[1];
    float  b  = *(float *)args[2];
    float  F  = *(float *)args[3];

    /* The output is a pointer to a complex float.                            */
    complex float *out = (complex float *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Float(x[i], a, b, F);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      complex_double_square_well                                            *
 *  Purpose:                                                                  *
 *      Same as complex_float_square_well but for doubles instead of floats.  *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) See complex_float_square_well above for more documentation.       *
 ******************************************************************************/
static void complex_double_square_well(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    double *x  =  (double *)args[0];
    double  a  = *(double *)args[1];
    double  b  = *(double *)args[2];
    double  F  = *(double *)args[3];

    /* The output is a pointer to a complex double.                           */
    complex double *out = (complex double *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      complex_long_double_square_well                                       *
 *  Purpose:                                                                  *
 *      Same as complex_float_square_well but for long doubles.               *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) See complex_float_square_well above for more documentation.       *
 ******************************************************************************/
static void complex_long_double_square_well(char **args, npy_intp *dimensions,
                                            npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    long double *x  =  (long double *)args[0];
    long double  a  = *(long double *)args[1];
    long double  b  = *(long double *)args[2];
    long double  F  = *(long double *)args[3];

    /* The output is a pointer to a complex double.                           */
    complex long double *out = (complex long double *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Long_Double(x[i], a, b, F);
    }
}

/*--------Inverted Square Well Diffraction Using Fresnel Approximation--------*/

static void complex_float_inv_square_well(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];
    float *x  =  (float *)args[0];
    float a   = *(float *)args[1];
    float b   = *(float *)args[2];
    float F   = *(float *)args[3];

    complex float *out = (complex float *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Float(x[i], a, b, F);
    }
}

static void complex_double_inv_square_well(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    double *x  =  (double *)args[0];
    double a   = *(double *)args[1];
    double b   = *(double *)args[2];
    double F   = *(double *)args[3];
    complex double *out = (complex double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

static void complex_long_double_inv_square_well(char **args,
                                                npy_intp *dimensions,
                                                npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    long double *x  =  (long double *)args[0];
    long double a   = *(long double *)args[1];
    long double b   = *(long double *)args[2];
    long double F   = *(long double *)args[3];
    complex long double *out = (complex long double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

/*-------------Phase from Square Well Using Fresnel Approximation-------------*/

static void float_square_well_phase(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    float *x   =  (float *)args[0];
    float a    = *(float *)args[1];
    float b    = *(float *)args[2];
    float F    = *(float *)args[3];
    float *out =  (float *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Float(x[i], a, b, F);
    }
}

static void double_square_well_phase(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    double *x   =  (double *)args[0];
    double a    = *(double *)args[1];
    double b    = *(double *)args[2];
    double F    = *(double *)args[3];
    double *out =  (double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Double(x[i], a, b, F);
    }
}

static void long_double_square_well_phase(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    long double *x   =  (long double *)args[0];
    long double a    = *(long double *)args[1];
    long double b    = *(long double *)args[2];
    long double F    = *(long double *)args[3];
    long double *out =  (long double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Double(x[i], a, b, F);
    }
}

PyUFuncGenericFunction sqwellsol_funcs[3] = {
    &complex_float_square_well,
    &complex_double_square_well,
    &complex_long_double_square_well
};

PyUFuncGenericFunction invsqwellsol_funcs[3] = {
    &complex_float_inv_square_well,
    &complex_double_inv_square_well,
    &complex_long_double_inv_square_well
};

PyUFuncGenericFunction sqwellphase_funcs[3] = {
    &float_square_well_phase,
    &double_square_well_phase,
    &long_double_square_well_phase
};

static void *PyuFunc_data[3
] = {NULL, NULL, NULL};

/*  Input and return types for square_well_diffraction.                       */
static char square_well_double_types[15] = {NPY_FLOAT,
                                            NPY_FLOAT,
                                            NPY_FLOAT,
                                            NPY_FLOAT, NPY_CFLOAT,
                                            NPY_DOUBLE,
                                            NPY_DOUBLE,
                                            NPY_DOUBLE,
                                            NPY_DOUBLE, NPY_CDOUBLE,
                                            NPY_LONGDOUBLE,
                                            NPY_LONGDOUBLE,
                                            NPY_LONGDOUBLE,
                                            NPY_LONGDOUBLE, NPY_CLONGDOUBLE};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fresnel_diffraction",
    NULL,
    -1,
    _fresnel_diffraction_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__fresnel_diffraction(void)
{
    PyObject *square_well_diffraction;
    PyObject *inverse_square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *m, *d;

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_data, square_well_double_types,
        3, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    inverse_square_well_diffraction = PyUFunc_FromFuncAndData(
        invsqwellsol_funcs, PyuFunc_data, square_well_double_types,
        3, 4, 1, PyUFunc_None, "inverse_square_well_diffraction", 
        "inverse_square_well_diffraction_docstring", 0
    );

    square_well_phase = PyUFunc_FromFuncAndData(
        sqwellphase_funcs, PyuFunc_data, square_well_double_types,
        3, 4, 1, PyUFunc_None, "square_well_phase", 
        "square_well_phase_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);

    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(square_well_phase);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *square_well_diffraction;
    PyObject *inverse_square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _fresnel_diffraction_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

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

    square_well_phase = PyUFunc_FromFuncAndData(
        sqwellphase_funcs, PyuFunc_data, sqwellsol_types,
        1, 4, 1, PyUFunc_None, "square_well_phase", 
        "square_well_phase_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);

    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(square_well_phase);

    return m;
}
#endif