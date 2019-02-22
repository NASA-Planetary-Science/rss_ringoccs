#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include <string.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

double SQRT_PI_BY_8 = 0.626657068657750125603941;
double HALF_PI = 1.5707963267948966;
double ONE_PI = 3.141592653589793;

static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

static void __rect(char *wfunc, double w_width, double dx,
                   npy_intp w_steps, char *nw_pts)
{
        int i;

        /* Window functions have an odd number of points. */
        *((int *)nw_pts) = (2 * (int)(w_width / (2.0 * dx))) + 1;

        for (i=0; i<*(int *)nw_pts; i++){
            *((double *)(wfunc+i*w_steps)) = 1.0;
        }
}

static void __coss(char *wfunc, double w_width, double dx,
                   npy_intp w_steps, char *nw_pts)
{
        int i;
        double x;

        /* Window functions have an odd number of points. */
        *(int *)nw_pts = (2 * (int)(w_width / (2.0 * dx))) + 1;

        for (i=0; i<*(int *)nw_pts; i++){
            x = i - (*(int *)nw_pts - 1) / 2.0;
            x *= ONE_PI * dx / w_width;
            x = cos(x);
            *((double *)(wfunc+i*w_steps)) = x*x;
        }
}

static void get_arr(char *x_arr, char *rho,
                      char *nw_pts, npy_intp rho_steps)
{
    int i, j;
    double x;
    i = -(*(int *)nw_pts-1)/2 * rho_steps;

    for (j=0; j<*(int *)nw_pts; j++){
        x = *(double *)(rho) - *(double *)(rho - i);
        x *= HALF_PI*x;
        *((double *)(x_arr+j*rho_steps)) = x;
        i += rho_steps;
    }
}

static void Fresnel_Transform_Func(char **args, npy_intp *dimensions,
                                      npy_intp* steps, void* data)
{
    char *n_pts = 0;
    double w_init, x, F, F2;
    static void (*fw)(char *, double, double, npy_intp, char *);
    npy_intp i, j, k, n_total;

    char *T_in              = args[0];
    char *rho_km_vals       = args[1];
    char *F_km_vals         = args[2];
    char *w_km_vals         = args[3];
    char *w_func            = args[4];
    char *x_arr             = args[5];
    char *start             = args[6];
    char *n_used            = args[7];
    char *dx_km             = args[8];
    char *wtype             = args[9];
    char *T_out             = args[10];

    npy_intp T_in_steps     = steps[0];
    npy_intp rho_steps      = steps[1];
    npy_intp F_steps        = steps[2];
    npy_intp w_steps        = steps[3];
    npy_intp w_f_steps      = steps[4];
    npy_intp x_arr_steps    = steps[5];
    npy_intp T_out_steps    = steps[9];

    n_total = *(int *)n_used * rho_steps;

    if (strncmp(*(char *)wtype, "rect", 4) == 0){
        fw = &__rect;
    }
    else if (strncmp(*(char *)wtype, "coss", 4) == 0){
        fw = &__coss;
    }
    else if (strncmp(*(char *)wtype, "kb20", 4) == 0){
        fw = &__coss;
    }
    else if (strncmp(*(char *)wtype, "kb25", 4) == 0){
        fw = &__coss;
    }
    else if (strncmp(*(char *)wtype, "kb35", 4) == 0){
        fw = &__coss;
    }
    else if (strncmp(*(char *)wtype, "kbmd20", 4) == 0){
        fw = &__coss;
    }
    else {
        fw = &__coss;
    }

    /* Compute first window width and window function. */
    w_km_vals   += *(int *)start*w_steps;
    rho_km_vals += *(int *)start*rho_steps;
    F_km_vals   += *(int *)start*F_steps;
    T_in        += *(int *)start*T_in_steps;
    T_out       += *(int *)start*T_out_steps;
    w_init       = *(double *)w_km_vals;

    fw(&w_func, w_init, *(double *)dx_km, w_f_steps, &n_pts);
    get_arr(&x_arr, &rho_km_vals, &n_pts, rho_steps);

    for (i=0; i<=n_total; ++i){
        /* Rho, Window width, Frensel scale for current point. */
        F = *(double *)F_km_vals;
        F2 = F * F;

        if (fabs(w_init - *(double *)w_km_vals) >= 2.0 * *(double *)dx_km) {
            /* Reset w_init and recompute window function. */
            w_init = *(double *)w_km_vals;
            fw(*w_func, w_init, *(double *)dx_km, w_f_steps, *n_pts);

            /* Compute psi for with stationary phase value */
            get_arr(x_arr, rho_km_vals, n_pts, rho_steps);
        }

        *((complex double *)T_out) = 0.0;
        k = - (*(int *)n_pts - 1) / 2;
        for (j=0; j<*(int *)n_pts; ++j){
            x = *(double *)(x_arr + j*x_arr_steps)/F2;
            *((complex double *)T_out) += *(double *)(w_func +j*w_f_steps) * (
                cos(x) - sin(x) * _Complex_I
            )* *(complex double *)(T_in + k*T_in_steps);
            k += 1;
        }

        *(complex double *)T_out *= *(double *)dx_km*(0.5+0.5*_Complex_I)/F;

        rho_km_vals += rho_steps;
        w_km_vals   += w_steps;
        F_km_vals   += F_steps;
        T_in        += T_in_steps;
        T_out       += T_out_steps; 
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_transform_funcs[1] = {&Fresnel_Transform_Func};

/* Input and return types for double input and out.. */
static char data_types[10] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_DOUBLE,
    NPY_COMPLEX128
};

static void *PyuFunc_data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_diffraction_functions",
    NULL,
    -1,
    _diffraction_functions_methods,
    NULL,
    NULL,
    NULL,
    NULL
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

    fresnel_transform = PyUFunc_FromFuncAndData(fresnel_transform_funcs, 
                                                PyuFunc_data, data_types,
                                                1, 9, 1, PyUFunc_None,
                                                "fresnel_transform",
                                                "fresnel_transform_docstring",
                                                0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform", fresnel_transform);
    Py_DECREF(fresnel_transform);
    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_transform;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _diffraction_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_transform = PyUFunc_FromFuncAndData(fresnel_transform_funcs, 
                                                PyuFunc_data, data_types,
                                                1, 9, 1, PyUFunc_None,
                                                "fresnel_transform",
                                                "fresnel_transform_docstring",
                                                0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform", fresnel_transform);
    Py_DECREF(fresnel_transform);
    return m;
}
#endif