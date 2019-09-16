/******************************************************************************
 *                          Diffraction Functions                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains wrappers for the routines found in                 *
 *      __diffraction_functions.c to allow use with the Python interpreter.   *
 ******************************************************************************
 *                             A FRIENDY WARNING                              *
 ******************************************************************************
 *  1.) This code uses complex numbers throughout, and is compatible with the *
 *      C99 standard. To use this code, make sure your compiler supports C99  *
 *      or more recent standards of the C Programming Language.               *
 * 
 *  2.) This code acts as Python wrappers for the pure C functions found in   *
 *      __diffraction_correction.h. As such, there is usage of the C-Numpy    *
 *      UFuncs API, as well as the standard C-Python API. This allows numpy   *
 *      arrays and various Python objects to be passed into these C routines. *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  Various header files required for the C-Python API to work.               */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  Window functions and Fresnel transforms defined here.                     */
#include "__diffraction_functions.h"

static void complex_double_fresnel_transform(char **args, npy_intp *dimensions,
                                             npy_intp *steps, void *data)
{
    DLPObj dlp;

    dlp.T_in         =  (complex double *)args[0];
    dlp.rho_km_vals  =  (double *)args[1];
    dlp.F_km_vals    =  (double *)args[2];
    dlp.phi_rad_vals =  (double *)args[3];
    dlp.kd_vals      =  (double *)args[4];
    dlp.B_rad_vals   =  (double *)args[5];
    dlp.D_km_vals    =  (double *)args[6];
    dlp.w_km_vals    =  (double *)args[7];
    dlp.start        = *(long *)args[8];
    dlp.n_used       = *(long *)args[9];
    dlp.wtype        = *(unsigned char *)args[10];
    dlp.use_norm     = *(unsigned char *)args[11];
    dlp.use_fwd      = *(unsigned char *)args[12];
    dlp.order        = *(unsigned char *)args[13];
    dlp.ecc          = *(double *)args[14];
    dlp.peri         = *(double *)args[15];
    dlp.T_out        =  (complex double *)args[16];

    if (dlp.order == 0){
        if ((dlp.ecc == 0.0) && (dlp.peri == 0.0)){
            DiffractionCorrectionNewton(dlp);
        }
        else {
            DiffractionCorrectionEllipse(dlp);
        }
    }
    else if (dlp.order == 1){
        DiffractionCorrectionFresnel(dlp);
    }
    else {
        DiffractionCorrectionLegendre(dlp);
    }
}

static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

/* Define pointers to the C functions. */
PyUFuncGenericFunction funcs[1] = {&complex_double_fresnel_transform};

/* Input and return types for the Fresnel Transform                           */
static char data_types[17] = {
    NPY_CDOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_UBYTE,
    NPY_UBYTE,
    NPY_UBYTE,
    NPY_UBYTE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_CDOUBLE
};

static void *PyuFunc_data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_diffraction_functions", NULL, -1,
    _diffraction_functions_methods, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit__diffraction_functions(void)
{
    PyObject *fresnel_transform;
    PyObject *m, *d;

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_transform = PyUFunc_FromFuncAndData(
        funcs, PyuFunc_data, data_types, 1, 16, 1, PyUFunc_None,
        "fresnel_transform", "fresnel_transform_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform", fresnel_transform);
    Py_DECREF(fresnel_transform);
    return m;
}
#else
PyMODINIT_FUNC init__diffraction_functions(void)
{
    PyObject *fresnel_transform;
    PyObject *m, *d;

    m = Py_InitModule("_diffraction_functions", _diffraction_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_transform = PyUFunc_FromFuncAndData(
        funcs, PyuFunc_data, data_types, 1, 16, 1, PyUFunc_None,
        "fresnel_transform_", "fresnel_transform_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform", fresnel_transform);
    Py_DECREF(fresnel_transform);
    return m;
}
#endif