#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

/*  For parallel processing.                              */
#include <omp.h>

/*  Various coefficients and constants defined here.      */
#include "__math_constants.h"

/*  Window functions and Fresnel Transforms defined here. */
#include "__window_functions.h"
#include "__diffraction_functions.h"


static void Fresnel_Transform_Quadratic_Func(char **args, npy_intp *dimensions,
                                             npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double,
                            double, long, npy_intp);

    char *T_in              = args[0];
    char *rho_km_vals       = args[1];
    char *F_km_vals         = args[2];
    char *w_km_vals         = args[3];
    char *start             = args[4];
    char *n_used            = args[5];
    char *wtype             = args[6];
    char *use_norm          = args[7];
    char *T_out             = args[8];

    npy_intp T_in_steps     = steps[0];
    npy_intp rho_steps      = steps[1];
    npy_intp F_steps        = steps[2];
    npy_intp w_steps        = steps[3];
    npy_intp T_out_steps    = steps[8];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else {fw = &Modified_Kaiser_Bessel_Window_2_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_transform;}
    else {FresT = &_fresnel_transform_norm;}

    /* Compute first window width and window function. */
    w_km_vals   += *(long *)start*w_steps;
    rho_km_vals += *(long *)start*rho_steps;
    F_km_vals   += *(long *)start*F_steps;
    T_in        += *(long *)start*T_in_steps;
    T_out       += *(long *)start*T_out_steps;
    w_init       = *(double *)w_km_vals;

    dx      = *(double *)(rho_km_vals+rho_steps) - *(double *)(rho_km_vals);
    two_dx  = 2.0*dx;
    nw_pts  = (int)(w_init / two_dx);

    double* x_arr   = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func  = (double *)malloc(sizeof(double)*nw_pts);

    get_arr(x_arr, dx, nw_pts);

    /* Compute Window Functions, and compute pi/2 * x^2   */
    for(j=0; j<nw_pts; ++j){
        w_func[j] = fw(x_arr[j], w_init);
        x_arr[j] *= PI_BY_TWO*x_arr[j];
    }

    #pragma omp parallel for
    for (i=0; i<=*(long *)n_used; ++i){
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            /* Reset w_init and recompute window function.    */
            w_init  = *(double *)w_km_vals;
            nw_pts  = (int)(w_init / two_dx);
            w_func  = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr   = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            get_arr(x_arr, dx, nw_pts);
            for(j=0; j<nw_pts; ++j){
                w_func[j] = fw(x_arr[j], w_init);
                x_arr[j] *= PI_BY_TWO*x_arr[j];
            }
        }

        *((complex double *)T_out) = FresT(x_arr, T_in, w_func,
                                           *(double *)F_km_vals, dx, nw_pts,
                                           T_in_steps);

        w_km_vals   += w_steps;
        F_km_vals   += F_steps;
        T_in        += T_in_steps;
        T_out       += T_out_steps;
    }
}

static void Fresnel_Transform_Quartic_Func(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, P_3, b_0, b_1, b_2;
    double A_0, A_1, A_2, rcpr_D, rcpr_F;
    double (*fw)(double, double);

    char *T_in              = args[0];
    char *F_km_vals         = args[2];
    char *phi_rad_vals      = args[3];
    char *kd_vals           = args[4];
    char *B_rad_vals        = args[5];
    char *D_km_vals         = args[6];
    char *w_km_vals         = args[7];
    char *start             = args[8];
    char *n_used            = args[9];
    char *wtype             = args[10];
    char *T_out             = args[11];

    npy_intp T_in_steps     = steps[0];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[11];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else {fw = &Modified_Kaiser_Bessel_Window_2_5;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    w_init  = *(double *)w_km_vals;
    dx      = *(double *)args[1];
    two_dx  = 2.0*dx;
    nw_pts  = (int)(w_init / (2.0 * dx));

    double* x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double* w_func  = (double *)malloc(sizeof(double) * nw_pts);

    get_arr(x_arr, dx, nw_pts);
    for(j=0; j<nw_pts; ++j){
        w_func[j] = fw(x_arr[j], w_init);
    }

    #pragma omp parallel for
    for (i=0; i<=*(long *)n_used; ++i){
        rcpr_F          = 1.0 / *(double *)F_km_vals;
        rcpr_D          = 1.0 / *(double *)D_km_vals;
        cosb            = cos(*(double *)B_rad_vals);
        cosp            = cos(*(double *)phi_rad_vals);
        sinp            = sin(*(double *)phi_rad_vals);
        Legendre_Coeff  = cosb*sinp;
        Legendre_Coeff *= Legendre_Coeff;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        /* Compute Legendre Polynomials. Rather than computing in the         *
         * standard form, i.e. (3 x^2 - 1) / 2, pre-compute the coefficients  *
         * of the various divisions, i.e. 1.5 x^2 - 0.5.                      */
        P_1 = cosb*cosp;
        P12 = P_1*P_1;
        P_2 = 1.5*P12-0.5;
        P_3 = (2.5*P12-1.5)*P_1;

        /* Second set of polynomials, (P_{n} - P_1 * P_{n+1}) / (n+2) */
        b_0 = 0.5 - 0.5*P12;
        b_1 = 0.333333333333 - 0.333333333333*P_2;
        b_2 = (P_2-P_1*P_3)*0.25;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = P_1*(b_1 - Legendre_Coeff*2.0*P_2);
        A_2 = b_2 - Legendre_Coeff*P_2*P_2;

        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = *(double *)w_km_vals;
            nw_pts  = (int)(w_init / two_dx);
            w_func  = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr   = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            get_arr(x_arr, dx, nw_pts);
            for(j=0; j<nw_pts; ++j){
                w_func[j] = fw(x_arr[j], w_init);
            }
        }

        /*  Compute the fresnel tranform about the current point.   */
        *((complex double *)T_out) = _fresnel_quartic(x_arr, T_in, w_func,
                                                      rcpr_D, A_0, A_1, A_2, dx,
                                                      rcpr_F, *(double *)kd_vals,
                                                      nw_pts, T_in_steps);
        
        /*  Increment pointers using pointer arithmetic, equivalent to        *
         *  changing var[n] to var[n+1].                                      */
        phi_rad_vals    += phi_steps;
        kd_vals         += kd_steps;
        B_rad_vals      += B_steps;
        D_km_vals       += D_steps;
        F_km_vals       += F_steps;
        w_km_vals       += w_steps;
        T_in            += T_in_steps;
        T_out           += T_out_steps;
    }
}

/*                            C-Python API Stuff                              */
static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

/* Define pointers to the C functions. */
PyUFuncGenericFunction f_quad_funcs[1] = {&Fresnel_Transform_Quadratic_Func};
PyUFuncGenericFunction f_qurt_funcs[1] = {&Fresnel_Transform_Quartic_Func};

/* Input and return types for Quadratic Fresnel Transform */
static char quad_data_types[9] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_COMPLEX128
};

/* Input and return types for Quartic Fresnel Transform */
static char qurt_data_types[12] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
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
static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT,
                                       "_diffraction_functions",
                                       NULL,
                                       -1,
                                       _diffraction_functions_methods,
                                       NULL, NULL, NULL, NULL};

PyMODINIT_FUNC PyInit__diffraction_functions(void)
{
    PyObject *fresnel_transform_quadratic;
    PyObject *fresnel_transform_quartic;

    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_transform_quadratic = PyUFunc_FromFuncAndData(
        f_quad_funcs, PyuFunc_data, quad_data_types,
        1, 8, 1, PyUFunc_None, "fresnel_transform_quadratic",
        "fresnel_transform_quadratic_docstring", 0
    );

    fresnel_transform_quartic = PyUFunc_FromFuncAndData(
        f_qurt_funcs, PyuFunc_data, qurt_data_types,
        1, 11, 1, PyUFunc_None, "fresnel_transform_quartic",
        "fresnel_transform_quartic_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform_quadratic",
                         fresnel_transform_quadratic);
    PyDict_SetItemString(d, "fresnel_transform_quartic",
                         fresnel_transform_quartic);
    Py_DECREF(fresnel_transform_quadratic);
    Py_DECREF(fresnel_transform_quartic);
    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_transform_quadratic;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _diffraction_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_transform_quadratic = PyUFunc_FromFuncAndData(
        f_quad_funcs, PyuFunc_data, quad_data_types,
        1, 8, 1, PyUFunc_None, "fresnel_transform_quadratic",
        "fresnel_transform_quadratic_docstring", 0
    );

    fresnel_transform_quartic = PyUFunc_FromFuncAndData(
        f_qurt_funcs, PyuFunc_data, qurt_data_types,
        1, 9, 1, PyUFunc_None, "fresnel_transform_quartic",
        "fresnel_transform_quartic_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform_quadratic",
                         fresnel_transform_quadratic);
    PyDict_SetItemString(d, "fresnel_transform_quartic",
                         fresnel_transform_quartic);
    Py_DECREF(fresnel_transform_quadratic);
    Py_DECREF(fresnel_transform_quartic);
    return m;
}
#endif