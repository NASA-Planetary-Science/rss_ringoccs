#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

/* Various constants used for diffraction correction. */
#define SQRT_PI_BY_8 0.626657068657750125603941
#define HALF_PI 1.5707963267948966
#define ONE_PI 3.141592653589793

/* Taylor Expansion of the Modified Bessel Function of the First Kind. */
#define MODIFIED_KAISER_BESSEL_2_5_A00 0.0
#define MODIFIED_KAISER_BESSEL_2_5_A01 0.04145269257683473
#define MODIFIED_KAISER_BESSEL_2_5_A02 0.1598131551145829
#define MODIFIED_KAISER_BESSEL_2_5_A03 0.27383552414424717
#define MODIFIED_KAISER_BESSEL_2_5_A04 0.26393049748717506
#define MODIFIED_KAISER_BESSEL_2_5_A05 0.16280559997382046
#define MODIFIED_KAISER_BESSEL_2_5_A06 0.06974074939338591
#define MODIFIED_KAISER_BESSEL_2_5_A07 0.02194877573814514
#define MODIFIED_KAISER_BESSEL_2_5_A08 0.005288714199797866
#define MODIFIED_KAISER_BESSEL_2_5_A09 0.0010068965459573457
#define MODIFIED_KAISER_BESSEL_2_5_A10 0.00015527610283472333
#define MODIFIED_KAISER_BESSEL_2_5_A11 1.978969147626019e-05
#define MODIFIED_KAISER_BESSEL_2_5_A12 2.1193188594867417e-06
#define MODIFIED_KAISER_BESSEL_2_5_A13 1.9338793216440484e-07
#define MODIFIED_KAISER_BESSEL_2_5_A14 1.521573968755883e-08
#define MODIFIED_KAISER_BESSEL_2_5_A15 1.0428703568483377e-09
#define MODIFIED_KAISER_BESSEL_2_5_A16 6.282176430490713e-11
#define MODIFIED_KAISER_BESSEL_2_5_A17 3.352216487175118e-12
#define MODIFIED_KAISER_BESSEL_2_5_A18 1.595536776390232e-13
#define MODIFIED_KAISER_BESSEL_2_5_A19 6.815840023528809e-15

static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

static void get_arr(double* x_arr, double dx, long nw_pts)
{
    /***************************************************************************
     *  Function:                                                              *
     *      get_arr                                                            *
     *  Purpose:                                                               *
     *      This computes the array pi/2 x^2, where x range from rho-w/2 to    *
     *      rho+w/2. Do the the symmetry involved in the computation of the    *
     *      Fresnel transform, on the points from rho-w/2 <= x < 0 are         *
     *      computed. The midpoint is always 0, and the right side is equal to *
     *      the left side. Hence, given a window with 2n+1 points, the array   *
     *      passed to get_arr will have n points. This cuts the number of      *
     *      computations in half while computing the same result.              *
     **************************************************************************/
    long i;
    double x;

    for (i=0; i<nw_pts; ++i){
        x = (i-nw_pts)*dx;
        x *= HALF_PI*x;
        x_arr[i] = x;
    }
}

static void get_arr_quartic(double* x_arr, double* x_arr_2,
                            double dx, long nw_pts)
{
    /***************************************************************************
     *  Function:                                                              *
     *      get_arr_quartic                                                    *
     *  Purpose:                                                               *
     *      This computes the array pi/2 x^2, where x range from rho-w/2 to    *
     *      rho+w/2. Do the the symmetry involved in the computation of the    *
     *      Fresnel transform, on the points from rho-w/2 <= x < 0 are         *
     *      computed. The midpoint is always 0, and the right side is equal to *
     *      the left side. Hence, given a window with 2n+1 points, the array   *
     *      passed to get_arr will have n points. This cuts the number of      *
     *      computations in half while computing the same result.              *
     **************************************************************************/
    long i;

    for (i=0; i<nw_pts; ++i){
        x_arr[i]    = (i-nw_pts)*dx;
        x_arr_2[i]  = x_arr[i]*x_arr[i];
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
        long i;
        double x;
        double bessel_x;
        dx = 2.0 * dx / w_width;

        for (i=0; i<nw_pts; i++){
            x = (i-nw_pts) * dx;
            x = 1.0 - x*x;
            bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
            bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;

            wfunc[i] = x*bessel_x;
        }
}

complex double _fresnel_transform(double* x_arr, char* T_in, double* w_func,
                                  double F, double dx, long n_pts,
                                  npy_intp T_in_steps)
{
    long i, j;
    double x, F2;
    complex double T_out, exp_negative_ix;

    T_out = 0.0;

    j = -n_pts;

    /*  Division is more expensive than multiplication, so store the
        reciprical of F as a variable and compute with that.          */
    F = 1.0/F;
    F2 = F*F;

    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*F2;
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        T_out += *(complex double *)T_in;
        T_out *= (0.5+0.5*_Complex_I)*dx*F;
    return T_out;
}

complex double _fresnel_quartic(double* x_arr, double* x_arr_2, char* T_in,
                                double* w_func, double rcpr_D, double rcpr_D2,
                                double A_0, double A_1, double A_2, double dx,
                                double rcpr_F, double kd, long n_pts,
                                npy_intp T_in_steps)
{
    long i, j;
    double x, x2, psi;
    complex double T_out, exp_negative_psi;

    T_out = 0.0;
    j = -n_pts;

    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x_arr_2[i]*rcpr_D2;

        /* Use Horner's Method to compute Psi. */
        psi = A_2*x+A_1;
        psi = psi*x+A_0;
        psi *= x2*kd;

        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi)) * w_func[i];
        T_out += exp_negative_psi * (*(complex double *)(T_in + j*T_in_steps) +
                                     *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        T_out += *(complex double *)T_in;
        T_out *= (0.5+0.5*_Complex_I)*dx*rcpr_F;
    return T_out;
}

static void Fresnel_Transform_Quadratic_Func(char **args, npy_intp *dimensions,
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
    nw_pts = (int)(w_init / (2.0 * dx));

    double* x_arr = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);

    
    fw(w_func, w_init, dx, nw_pts);
    get_arr(x_arr, dx, nw_pts);

    for (i=0; i<=*(long *)n_used; ++i){
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init = *(double *)w_km_vals;
            nw_pts = (int)(w_init / (2.0 * dx));
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

static void Fresnel_Transform_Quartic_Func(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data)
{
    long i, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, P_3, b_0, b_1, b_2;
    double A_0, A_1, A_2, rcpr_D, rcpr_D2, rcpr_F;
    void (*fw)(double*, double, double, long);

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
    else if (*(int *)wtype == 2){fw = &__coss;}
    else if (*(int *)wtype == 3){fw = &__coss;}
    else if (*(int *)wtype == 4){fw = &__coss;}
    else if (*(int *)wtype == 5){fw = &__coss;}
    else {fw = &__kbmd25;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;
    w_init           = *(double *)w_km_vals;

    dx      = *(double *)args[1];
    two_dx  = 2.0*dx;
    nw_pts  = (int)(w_init / (2.0 * dx));

    double* x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double* x_arr_2 = (double *)malloc(sizeof(double) * nw_pts);
    double* w_func  = (double *)malloc(sizeof(double) * nw_pts);

    fw(w_func, w_init, dx, nw_pts);
    get_arr_quartic(x_arr, x_arr_2, dx, nw_pts);

    for (i=0; i<=*(long *)n_used; ++i){
        rcpr_F          = 1.0 / *(double *)F_km_vals;
        rcpr_D          = 1.0 / *(double *)D_km_vals;
        rcpr_D2         = rcpr_D*rcpr_D;
        cosb            = cos(*(double *)B_rad_vals);
        cosp            = cos(*(double *)phi_rad_vals);
        sinp            = sin(*(double *)phi_rad_vals);
        Legendre_Coeff  = cosb*cosb*sinp*sinp;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        /* Compute Legendre Polynomials. */
        P_1 = cosb*cosp;
        P12 = P_1*P_1;
        P_2 = (3.0*P12-1.0)*0.5;
        P_3 = (5.0*P12-3.0)*0.5*P_1;

        /* Second set of polynomials, (P_{n} - P_1 * P_{n+1}) / (n+2) */
        b_0 = 0.5 - 0.5*P12;
        b_1 = (P_1-P_1*P_2)*0.333333333333;
        b_2 = (P_2-P_1*P_3)*0.25;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = b_1 - Legendre_Coeff*2.0*P_1*P_2;
        A_2 = b_2 - Legendre_Coeff*P_2*P_2;

        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = *(double *)w_km_vals;
            nw_pts  = (int)(w_init / (2.0 * dx));
            w_func  = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr   = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            x_arr_2 = (double *)realloc(x_arr_2, sizeof(double)*nw_pts);
            fw(w_func, w_init, dx, nw_pts);
            get_arr_quartic(x_arr, x_arr_2, dx, nw_pts);
        }

        *((complex double *)T_out) = _fresnel_quartic(x_arr, x_arr_2, T_in,
                                                      w_func, rcpr_D, rcpr_D2,
                                                      A_0, A_1, A_2, dx, rcpr_F,
                                                      *(double *)kd_vals,
                                                      nw_pts, T_in_steps);
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

/* Define pointers to the C functions. */
PyUFuncGenericFunction f_quad_funcs[1] = {&Fresnel_Transform_Quadratic_Func};
PyUFuncGenericFunction f_qurt_funcs[1] = {&Fresnel_Transform_Quartic_Func};

/* Input and return types for Quadratic Fresnel Transform */
static char quad_data_types[8] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_COMPLEX128
};

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
        1, 7, 1, PyUFunc_None, "fresnel_transform_quadratic",
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
        1, 7, 1, PyUFunc_None, "fresnel_transform_quadratic",
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