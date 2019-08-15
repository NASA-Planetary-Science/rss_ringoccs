/*  To avoid compiler warnings about deprecated numpy stuff.                 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* cosine and sine are defined here. */
#include <math.h>

/*  complex data types, as well as _Complex_I, are defined here.             */
#include <complex.h>

/* Include fresnel integrals header. This includes frensel_sin/cos.          */
#include "__fresnel_integrals.h"
#include "__fresnel_kernel.h"

/*  Where compute_norm_eq lives.                                             */
#include "__normalized_equivalent_width.h"

/*  Various header files required for the C-Python API to work.              */
#include "../../../include/Python.h"
#include "../../../include/ndarraytypes.h"
#include "../../../include/ufuncobject.h"

/*---------------------------DEFINE PYTHON FUNCTIONS-------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above *
 * functions to be called in Python. Numpy arrays, as well as floating point *
 * and integer valued arguments may then be passed into these functions for  *
 * improvement in performance, as opposed to the routines written purely in  *
 * Python. Successful compiling requires the Numpy and Python header files.  *
 *---------------------------------------------------------------------------*/
static PyObject *compute_norm_eq(PyObject *self, PyObject *args)
{
    #define FNAME "rss_ringoccs.diffrec.special_functions.compute_norm_eq\n"
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long normeq;
        PyArg_ParseTuple(args, "l", &normeq);
        return PyLong_FromLong(normeq);
    }
    else if (PyFloat_Check(tuple)){
        double normeq;
        PyArg_ParseTuple(args, "d", &normeq);
        return PyFloat_FromDouble(normeq);
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput must be a one-dimensional array.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);
        if (dim == 0){
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput array is empty.");
            return NULL;
        }

        if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(Normeq_Double((double *)data, dim));
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(Normeq_Long_Double((long double *)data, dim));
        }
        else if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(Normeq_Float((float *)data, dim));
        }
        else if (typenum == NPY_LONG){
            return PyFloat_FromDouble(Normeq_Long((long *)data, dim));
        }
        else {
            PyErr_Format(PyExc_TypeError, FNAME
                         "\rInput should be a numpy array of real numbers.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError, FNAME
                     "\rInput should be a numpy array of real numbers.");
        return NULL;
    }
}

static PyMethodDef _special_functions_methods[] =
{
    {"compute_norm_eq", compute_norm_eq,
     METH_VARARGS, "Compute the maximum of a numpy array."},
    {NULL, NULL, 0, NULL}
};
/*-------------------------DEFINE UNIVERSAL FUNCTIONS-------------------------*/

/*  Functions from __fresnel_integrals.h                                      */
static void double_fresnelsin(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n  = dimensions[0];
    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Fresnel_Sine_Taylor_to_Asymptotic_Func(in[i]);
    }
}

static void double_fresnelcos(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n  = dimensions[0];
    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Cosine_Taylor_to_Asymptotic_Func is defined  *
         *  in _fresnel_cosine.h. Make sure this is in the current directory! */
        out[i] = Fresnel_Cosine_Taylor_to_Asymptotic_Func(in[i]);
    }
}

static void complex_sqwellsol(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    double *x =  (double *)args[0];
    double  a  = *(double *)args[1];
    double  b  = *(double *)args[2];
    double  F  = *(double *)args[3];
    double complex *out = (double complex *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Func(x[i], a, b, F);
    }
}

static void complex_invsqwellsol(char **args, npy_intp *dimensions,
                                 npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    double *x  =  (double *)args[0];
    double a   = *(double *)args[1];
    double b   = *(double *)args[2];
    double F   = *(double *)args[3];
    double complex *out = (double complex *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Func(x[i], a, b, F);
    }
}

static void double_sqwellphase(char **args, npy_intp *dimensions,
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
        *((double complex*)out) = Square_Well_Diffraction_Phase_Func(
            *(double *)x, *(double *)a, *(double *)b, *(double *)F
        );
        /*END main ufunc computation*/

        x += in_step;
        out += out_step;
    }
}

/*  Functions from __fresnel_kernel.h                                         */
static void double_psi(char **args, npy_intp *dimensions,
                       npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((double *)out) = Fresnel_Psi_Func(*(double *)kD, *(double *)rho,
                                            *(double *)rho0, *(double *)phi,
                                            *(double *)phi0, *(double *)B,
                                            *(double *)D);

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void double_dpsi_dphi(char **args, npy_intp *dimensions,
                             npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((double *)out) = Fresnel_dPsi_dPhi_Func(
            *(double *)kD, *(double *)rho, *(double *)rho0, *(double *)phi,
            *(double *)phi0, *(double *)B, *(double *)D
        );

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

/*  Define pointers to the C functions.                                       */
PyUFuncGenericFunction fresnel_sin_funcs[1]     = {&double_fresnelsin};
PyUFuncGenericFunction fresnel_cos_funcs[1]     = {&double_fresnelcos};
PyUFuncGenericFunction sqwellsol_funcs[1]       = {&complex_sqwellsol};
PyUFuncGenericFunction invsqwellsol_funcs[1]    = {&complex_invsqwellsol};
PyUFuncGenericFunction sqwellphase_funcs[1]     = {&double_sqwellphase};
PyUFuncGenericFunction psi_funcs[1]             = {&double_psi};
PyUFuncGenericFunction dpsi_funcs[1]            = {&double_dpsi_dphi};

/*  Input and return types for double input and out.                          */
static char double_double_types[2] = {NPY_DOUBLE, NPY_DOUBLE};
static void *PyuFunc_data[1] = {NULL};

/*  Input and return types for square_well_diffraction.                       */
static char sqwellsol_types[5] = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                  NPY_DOUBLE, NPY_COMPLEX128};

/* Input and return types for fresnel_psi.                                    */
static char octo_double_types[8] = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                    NPY_DOUBLE, NPY_DOUBLE};

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
    PyObject *inverse_square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *fresnel_psi;
    PyObject *fresnel_dpsi_dphi;
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

    fresnel_psi = PyUFunc_FromFuncAndData(
        psi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
        PyUFunc_None, "fresnel_psi",  "fresnel_psi_docstring", 0
    );

    fresnel_dpsi_dphi = PyUFunc_FromFuncAndData(
        dpsi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
        PyUFunc_None, "fresnel_dpsi_dphi",  "fresnel_dpsi_dphi_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);

    Py_DECREF(fresnel_sin);
    Py_DECREF(fresnel_cos);
    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(square_well_phase);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_sin;
    PyObject *fresnel_cos;
    PyObject *square_well_diffraction;
    PyObject *inverse_square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *fresnel_psi;
    PyObject *fresnel_dpsi_dphi;
    PyObject *m, *d;

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

    fresnel_psi = PyUFunc_FromFuncAndData(
        psi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
        PyUFunc_None, "fresnel_psi",  "fresnel_psi_docstring", 0
    );

    fresnel_dpsi_dphi = PyUFunc_FromFuncAndData(
        dpsi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
        PyUFunc_None, "fresnel_dpsi_dphi",  "fresnel_dpsi_dphi_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);

    Py_DECREF(fresnel_sin);
    Py_DECREF(fresnel_cos);
    Py_DECREF(square_well_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(square_well_phase);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);

    return m;
}
#endif