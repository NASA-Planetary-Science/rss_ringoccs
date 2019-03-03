#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

#define SQRT_PI_BY_8 0.626657068657750125603941
#define HALF_PI 1.5707963267948966
#define ONE_PI 3.141592653589793

/*
#define KAISER_BESSEL_2_5_A00 1.0
#define KAISER_BESSEL_2_5_A01 -14.71036758931379
#define KAISER_BESSEL_2_5_A02 94.28094137560409
#define KAISER_BESSEL_2_5_A03 -353.52377591977927
#define KAISER_BESSEL_2_5_A04 878.0028262813076
#define KAISER_BESSEL_2_5_A05 -1551.8157242842929
#define KAISER_BESSEL_2_5_A06 2048.5647177480982
*/

#define MODIFIED_KAISER_BESSEL_2_5_A00 0.0026880229613119944
#define MODIFIED_KAISER_BESSEL_2_5_A01 0.0414526925768658
#define MODIFIED_KAISER_BESSEL_2_5_A02 0.15981315511470268
#define MODIFIED_KAISER_BESSEL_2_5_A03 0.27383552414445234
#define MODIFIED_KAISER_BESSEL_2_5_A04 0.26393049748737285
#define MODIFIED_KAISER_BESSEL_2_5_A05 0.16280559997394245
#define MODIFIED_KAISER_BESSEL_2_5_A06 0.06974074939343818
#define MODIFIED_KAISER_BESSEL_2_5_A07 0.021948775738161588
#define MODIFIED_KAISER_BESSEL_2_5_A08 0.005288714199801829
#define MODIFIED_KAISER_BESSEL_2_5_A09 0.0010068965459581003
#define MODIFIED_KAISER_BESSEL_2_5_A10 0.00015527610283483967
#define MODIFIED_KAISER_BESSEL_2_5_A11 1.978969147627502e-05


static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

static void get_arr(double* x_arr, double dx, long nw_pts)
{
    long i, j;
    double x;

    j = -(nw_pts-1)/2;

    for (i=0; i<nw_pts; ++i){
        x = - j*dx;
        x *= HALF_PI*x;
        x_arr[i] = x;
        j += 1;
    }
}

static void __rect(double* wfunc, double w_width, double dx, long nw_pts)
{
        long i;

        for (i=0; i<nw_pts; i++){
            wfunc[i] = 1.0;
        }
}

static void __coss(double* wfunc, double w_width, double dx, long nw_pts)
{
        long i, j;
        double x;
        dx = ONE_PI * dx / w_width;

        j = -(nw_pts - 1) / 2.0;
        for (i=0; i<nw_pts; i++){
            x = j * dx;
            x = cos(x);
            x *= x;
            wfunc[i] = x;
            j += 1;
        }
}

static void __kbmd25(double* wfunc, double w_width, double dx, long nw_pts)
{
        long i, j;
        double x;
        double bessel_x;
        dx = 2.0 * dx / w_width;

        j = -(nw_pts - 1) / 2.0;
        for (i=0; i<nw_pts; i++){
            x = j * dx;
            x *= x;
            x = 1.0 - x;
            bessel_x = MODIFIED_KAISER_BESSEL_2_5_A08*x;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;

            wfunc[i] = x*bessel_x;
            j += 1;
        }
}

complex double _fresnel_transform(double* x_arr, char* T_in, double* w_func,
                                  double F, double dx, long n_pts,
                                  npy_intp T_in_steps)
{
    long i, j;
    double x, F2;
    complex double T_out = 0.0;

    j = -(n_pts - 1)/2;
    F = 1.0/F;
    F2 = F*F;

    for (i=0; i<n_pts; ++i){
        x = x_arr[i]*F2;
        T_out += w_func[i] * (cos(x) - sin(x) * _Complex_I) *
                 *(complex double *)(T_in + j*T_in_steps);
        j += 1;
    }

        T_out *= (0.5+0.5*_Complex_I)*dx*F;
    return T_out;
}

static void Fresnel_Transform_Func(char **args, npy_intp *dimensions,
                                   npy_intp* steps, void* data)
{
    long i, nw_pts;
    double w_init, dx, two_dx;
    void (*fw)(double*, double, double, long);

    char *T_in              = args[0];
    char *rho_km_vals       = args[1];
    char *F_km_vals         = args[2];
    char *w_km_vals         = args[3];
    char *start             = args[4];
    char *n_used            = args[5];
    char *wtype             = args[6];
    char *T_out             = args[7];

    npy_intp T_in_steps     = steps[0];
    npy_intp rho_steps      = steps[1];
    npy_intp F_steps        = steps[2];
    npy_intp w_steps        = steps[3];
    npy_intp T_out_steps    = steps[7];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &__coss;}
    else if (*(int *)wtype == 3){fw = &__coss;}
    else if (*(int *)wtype == 4){fw = &__coss;}
    else if (*(int *)wtype == 5){fw = &__coss;}
    else {fw = &__kbmd25;}

    /* Compute first window width and window function. */
    w_km_vals   += *(long *)start*w_steps;
    rho_km_vals += *(long *)start*rho_steps;
    F_km_vals   += *(long *)start*F_steps;
    T_in        += *(long *)start*T_in_steps;
    T_out       += *(long *)start*T_out_steps;
    w_init       = *(double *)w_km_vals;


    dx = *(double *)(rho_km_vals+rho_steps) - *(double *)(rho_km_vals);
    two_dx = 2.0*dx;
    nw_pts = (2 * (int)(w_init / (2.0 * dx))) + 1;

    double* x_arr = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);

    
    fw(w_func, w_init, dx, nw_pts);
    get_arr(x_arr, dx, nw_pts);

    for (i=0; i<=*(long *)n_used; ++i){
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init = *(double *)w_km_vals;
            nw_pts = (2 * (int)(w_init / (2.0 * dx))) + 1;
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            fw(w_func, w_init, dx, nw_pts);
            get_arr(x_arr, dx, nw_pts);
        }

        *((complex double *)T_out) = _fresnel_transform(x_arr, T_in, w_func,
                                                        *(double *)F_km_vals,
                                                        dx,  nw_pts,
                                                        T_in_steps);

        w_km_vals   += w_steps;
        F_km_vals   += F_steps;
        T_in        += T_in_steps;
        T_out       += T_out_steps;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_transform_funcs[1] = {&Fresnel_Transform_Func};

/* Input and return types for double input and out.. */
static char data_types[8] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
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
                                                1, 7, 1, PyUFunc_None,
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