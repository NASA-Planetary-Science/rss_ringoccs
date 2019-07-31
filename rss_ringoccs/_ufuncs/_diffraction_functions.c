/*******************************************************************************
 *                          Diffraction Functions                              *
 *******************************************************************************
 *  Purpose:                                                                   *
 *      This file contains functions used for computing the Fresnel Inverse    *
 *      Transform on a set of diffraction limited data. There are several      *
 *      Methods of performing this:                                            *
 *          Fresnel Quadratic Approximations:                                  *
 *              Classic quadratic approximation that is used in Fourier Optics.*
 *          Legendre Cubic Expansion:                                          *
 *              Approximation of the Fresnel kernel by a quartic polynomial    *
 *              using the Legendre Polynomials.                                *
 *          Legendre Quartic Expansion:                                        *
 *              Approximation of the Fresnel kernel by a quartic polynomial    *
 *              using the Legendre Polynomials.                                *
 *          Legendre Sextic Expansion:                                         *
 *              Approximation of the Fresnel kernel by a sextic polynomial     *
 *              using the Legendre Polynomials.                                *
 *          Legendre Octic Expansion:                                          *
 *              Approximation of the Fresnel kernel by a octic polynomial      *
 *              using the Legendre Polynomials.                                *
 *******************************************************************************
 *  Variables:                                                                 *
 *     A_0         (Double):                                                   *
 *         The coefficient of the x^2 term in the expansion for psi.           *
 *     A_1         (Double):                                                   *
 *         The coefficient of the x^3 term in the expansion for psi.           *
 *     A_2         (Double):                                                   *
 *         The coefficient of the x^4 term in the expansion for psi.           *
 *     A_3         (Double):                                                   *
 *         The coefficient of the x^5 term in the expansion for psi.           *
 *     A_4         (Double):                                                   *
 *         The coefficient of the x^6 term in the expansion for psi.           *
 *     A_5         (Double):                                                   *
 *         The coefficient of the x^7 term in the expansion for psi.           *
 *     A_6         (Double):                                                   *
 *         The coefficient of the x^8 term in the expansion for psi.           *
 *     dx          (Double):                                                   *
 *         Spacing between two consecutive ring intercept points, km.          *
 *     kd          (Double):                                                   *
 *         The wavenumber, k, weighted by the spacecraft-ring distance, D.     *
 *     n_pts       (Long):                                                     *
 *         Half the number of points in the window, rounded down.              *
 *     rcpr_D      (Double):                                                   *
 *         1/D, where D is the distant from the spacecraft to the ring         *
 *         intercept point, in kilometers.                                     *
 *     rcpr_F      (Double):                                                   *
 *         The reciprocal of the Fresnel scale in kilometers.                  *
 *     T_in        (Pointer to Char):                                          *
 *         The raw diffraction-limited data that is to be corrected.           *
 *     T_in_steps  (npy_intp (Equivalent to Long)):                            *
 *         The number of steps in memory to get from the nth data point        *
 *         to the (n+1)th data point in the T_in variable (See above).         *
 *     T_out       (Complex Double):                                           *
 *         The diffraction corrected profile.                                  *
 *     w_func      (Pointer to Double):                                        *
 *         Pre-computed window function. Should have the same number of        *
 *         points as the x_arr pointer.                                        *
 *     x_arr       (Pointer to Double):                                        *
 *         The ring radii within the given window.                             *
 *******************************************************************************
 *  The Inverse Fresnel Transform:                                             *
 *                                                                             *
 *                W/2                                                          *
 *                 -                                                           *
 *                | |                                                          *
 *     T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0                *
 *              | |                                                            *
 *               -                                                             *
 *              -W/2                                                           *
 *                                                                             *
 *  Where T_hat is the diffracted data, w is the window function, r is         *
 *  the ring intercept point, and r_0 is a dummy variable of integration.      *
 *  psi is the Fresnel Kernel, and exp is simply the exponential function.     *
 *******************************************************************************
 *  The Normalization Scheme:                                                  *
 *      As the resolution get's too high, say 10 km or greater, then window    *
 *      width quickly shrinks to zero. Thus the integral will be approximately *
 *      zero. To account for this, the option to normalize the integral by the *
 *      window width is offered. The normalization is defined as follows:      *
 *                                                                             *
 *                    |     _ +infinity           |                            *
 *                    |    | |                    |                            *
 *                    |    |    exp(-i psi(x)) dx |                            *
 *                    |  | |                      |                            *
 *                    |   -  -infinity            |                            *
 *          Norm =  __________________________________                         *
 *                  |    -  +W/2                    |                          *
 *                  |   | |                         |                          *
 *                  |   |    w(x) exp(-i psi(x)) dx |                          *
 *                  | | |                           |                          *
 *                  |  -   -W/2                     |                          *
 *                                                                             *
 *      This has the effect of making the free-space regions, or regions which *
 *      were not affected by diffraction, evaluate to approximately one,       *
 *      regardless of what resolution was chosen.                              *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  get_arr:                                                                   *
 *      Void function that takes in a pointer to a double array and creates    *
 *      an array of values for the ring radius within half a window width of   *
 *      the ring intercept point. Do to symmetry, only the values to the left  *
 *      of the ring intercept point are computed.                              *
 *******************************************************************************
 *  _fresnel_transform:                                                        *
 *      Computes the Fresnel Inverse Transform using the classic Fresnel       *
 *      quadratic approximation. No normalization is applied.                  *
 *******************************************************************************
 *  _fresnel_transform_norm:                                                   *
 *      Same as _fresnel_transform, but the normalization is applied.          *
 *******************************************************************************
 *  _fresnel_cubic,                                                            *
 *  _fresnel_quartic,                                                          *
 *  _fresnel_sextic,                                                           *
 *  _fresnel_octic:                                                            *
 *      Computes the Fresnel Inverse Transform using Legendre Polynomials to   *
 *      approximate the Fresnel kernel to various powers (3, 4, 6, or 8).      *
 *******************************************************************************
 *  _fresnel_cubic_norm,                                                       *
 *  _fresnel_quartic_norm,                                                     *
 *  _fresnel_sextic_norm,                                                      *
 *  _fresnel_octic_norm:                                                       *
 *      Same as previous functions, but with the normalization scheme.         *
 *******************************************************************************
 *                             A FRIENDY WARNING                               *
 *******************************************************************************
 *  This code uses complex numbers throughout, and is compatible with the C99  *
 *  standard. To use this code, make sure your compiler supports C99 or more   *
 *  recent standards of the C Programming Language.                            *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  Various trig functions, complex variables, and more found here.           */
#include <math.h>
#include <complex.h>

/*  Various header files required for the C-Python API to work.               */
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  Window functions and Fresnel Transforms defined here.                     */
#include "__window_functions.h"
#include "__diffraction_functions.h"

/*  Functions for computing the Fresnel Kernel and Newton's Method.           */
#include "__fresnel_kernel.h"


static void Fresnel_Transform_Quadratic_Func(char **args, npy_intp *dimensions,
                                             npy_intp* steps, void* data)
{
    /*  i and j used for indexing, nw_pts is number of points in window.      */
    long i, j, nw_pts;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).   */
    double w_init, dx, two_dx;

    /*  Declare function pointers. fw is window function, FresT is Fresnel    *
     *  transform. Options for FresT are with or without normalization.       */
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double,
                            double, long, npy_intp);

    /*  Grab all data passed in from Python.                                  */
    char *T_in           = args[0];
    char *F_km_vals      = args[2];
    char *w_km_vals      = args[3];
    char *start          = args[4];
    char *n_used         = args[5];
    char *wtype          = args[6];
    char *use_norm       = args[7];
    char *use_fwd        = args[8];
    char *T_out          = args[9];

    /*  Number of steps in memory from the ith point to the (i+1)th point.    */
    npy_intp T_in_steps  = steps[0];
    npy_intp F_steps     = steps[2];
    npy_intp w_steps     = steps[3];
    npy_intp T_out_steps = steps[9];

    /*  Cast the selected window type to the fw pointer.                      */
    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    /*  Cast FresT to the appropriate function.                               */
    if (*(int *)use_norm == 0){FresT = &_fresnel_transform;}
    else {FresT = &_fresnel_transform_norm;}

    /* Move the pointers to the correct starting point. Compute window width. */
    w_km_vals += *(long *)start*w_steps;
    F_km_vals += *(long *)start*F_steps;
    T_in      += *(long *)start*T_in_steps;
    T_out     += *(long *)start*T_out_steps;
    w_init     = *(double *)w_km_vals;

    /*  Compute some extra necessary variables.                               */
    dx     = *(double *)args[1];
    two_dx = 2.0*dx;
    nw_pts = (int)(w_init / (2.0 * dx));

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx.                                          */
    double* x_arr  = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);

    /*  Pass the x_arr array (ring radius) to the void function get_arr.      *
     *  This alters the x_arr pointer so that it's values range from -W/2 to  *
     *  zero, where W is the window width.                                    */
    get_arr(x_arr, dx, nw_pts);

    /* Compute Window Functions, and compute pi/2 * x^2                       */
    for(j=0; j<nw_pts; ++j){
        w_func[j] = fw(x_arr[j], w_init);
        x_arr[j] *= PI_BY_TWO*x_arr[j];

        /*  If forward transform is selected, negate x_arr.                   */
        x_arr[j] *= (*(int *)use_fwd == 0) - (*(int *)use_fwd == 1);
    }

    /*  Compute the Fresnel transform across the input data.                  */
    for (i=0; i<=*(long *)n_used; ++i){
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            /* Reset w_init and recompute window function.                    */
            w_init = *(double *)w_km_vals;
            nw_pts = (int)(w_init / two_dx);

            /*  Reallocate memory, since the sizes of the arrays changed.     */
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.         */
            get_arr(x_arr, dx, nw_pts);

            /* Compute Window Functions, and compute pi/2 * x^2               */
            for(j=0; j<nw_pts; ++j){
                w_func[j] = fw(x_arr[j], w_init);
                x_arr[j] *= PI_BY_TWO*x_arr[j];
                x_arr[j] *= (*(int *)use_fwd == 0) - (*(int *)use_fwd == 1);
            }
        }

        /*  Compute the Fresnel Transform about the current point.            */
        *((complex double *)T_out) = FresT(x_arr, T_in, w_func,
                                           *(double *)F_km_vals, dx, nw_pts,
                                           T_in_steps);

        /*  Move the pointers to the next point.                              */
        w_km_vals   += w_steps;
        F_km_vals   += F_steps;
        T_in        += T_in_steps;
        T_out       += T_out_steps;
    }
    free(x_arr);
    free(w_func);
}

static void Fresnel_Transform_Cubic_Func(char **args, npy_intp *dimensions,
                                         npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, b_0, b_1;
    double A_0, A_1, rcpr_D, rcpr_F;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double, double, double,
                            double, double, double, long, npy_intp);

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
    char *use_norm          = args[11];
    char *use_fwd           = args[12];
    char *T_out             = args[13];

    npy_intp T_in_steps     = steps[0];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[13];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_cubic;}
    else {FresT = &_fresnel_cubic_norm;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    if (*(int *)use_fwd == 1){
        for (i=0; i<=*(long *)n_used; ++i){
            *(double *)(kd_vals+i*kd_steps) *= -1.0;
        }
    }

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

        /* Second set of polynomials, (P_{n} - P_1 * P_{n+1}) / (n+2) */
        b_0 = 0.5 - 0.5*P12;
        b_1 = (0.333333333333 - 0.333333333333*P_2)*P_1;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = b_1 - Legendre_Coeff*2.0*P_1*P_2;

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
        *((complex double *)T_out) = FresT(x_arr, T_in, w_func, rcpr_D,
                                           A_0, A_1, dx, rcpr_F,
                                           *(double *)kd_vals,
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
    free(x_arr);
    free(w_func);
}

static void Fresnel_Transform_Quartic_Func(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, P_3, b_0, b_1, b_2;
    double A_0, A_1, A_2, rcpr_D, rcpr_F;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double, double, double,
                            double, double, double, double, long, npy_intp);

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
    char *use_norm          = args[11];
    char *use_fwd           = args[12];
    char *T_out             = args[13];

    npy_intp T_in_steps     = steps[0];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[13];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_quartic;}
    else {FresT = &_fresnel_quartic_norm;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    if (*(int *)use_fwd == 1){
        for (i=0; i<=*(long *)n_used; ++i){
            *(double *)(kd_vals+i*kd_steps) *= -1.0;
        }
    }

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
        b_1 = (0.333333333333 - 0.333333333333*P_2)*P_1;
        b_2 = (P_2 - P_1*P_3)*0.25;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = b_1 - Legendre_Coeff*2.0*P_1*P_2;
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
        *((complex double *)T_out) = FresT(x_arr, T_in, w_func, rcpr_D,
                                           A_0, A_1, A_2, dx, rcpr_F,
                                           *(double *)kd_vals,
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
    free(x_arr);
    free(w_func);
}

static void Fresnel_Transform_Sextic_Func(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, P_3, P_4, P_5;
    double b_0, b_1, b_2, b_3, b_4;
    double A_0, A_1, A_2, A_3, A_4, rcpr_D, rcpr_F;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double, double, double,
                            double, double, double, double, double, double,
                            long, npy_intp);

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
    char *use_norm          = args[11];
    char *use_fwd           = args[12];
    char *T_out             = args[13];

    npy_intp T_in_steps     = steps[0];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[13];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_sextic;}
    else {FresT = &_fresnel_sextic_norm;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    if (*(int *)use_fwd == 1){
        for (i=0; i<=*(long *)n_used; ++i){
            *(double *)(kd_vals+i*kd_steps) *= -1.0;
        }
    }

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
        P_4 = (4.375*P12-7.5)*P12+0.375;
        P_5 = ((7.875*P12-8.75)*P12+1.875)*P_1;

        /* Second set of polynomials, (P_{n} - P_1 * P_{n+1}) / (n+2) */
        b_0 = 0.5 - 0.5*P12;
        b_1 = (0.333333333333 - 0.333333333333*P_2)*P_1;
        b_2 = (P_2 - P_1*P_3)*0.25;
        b_3 = (P_3 - P_1*P_4)*0.20;
        b_4 = (P_4 - P_1*P_5)*0.16666666666666666;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = b_1 - Legendre_Coeff*2.0*P_1*P_2;
        A_2 = b_2 - Legendre_Coeff*P_2*P_2;
        A_3 = b_3 - Legendre_Coeff*2.0*(P_1*P_4+P_2*P_3);
        A_4 = b_4 - Legendre_Coeff*(2.0*P_2*P_4+P_3*P_3);

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
        *((complex double *)T_out) = FresT(x_arr, T_in, w_func, rcpr_D,
                                           A_0, A_1, A_2, A_3, A_4, dx, rcpr_F,
                                           *(double *)kd_vals,
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
    free(x_arr);
    free(w_func);
}

static void Fresnel_Transform_Octic_Func(char **args, npy_intp *dimensions,
                                         npy_intp* steps, void* data)
{
    long i, j, nw_pts;
    double w_init, dx, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, P_1, P12, P_2, P_3, P_4, P_5, P_6, P_7;
    double b_0, b_1, b_2, b_3, b_4, b_5, b_6;
    double A_0, A_1, A_2, A_3, A_4, A_5, A_6;
    double rcpr_D, rcpr_F;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double, double, double,
                            double, double, double, double, double, double,
                            double, double, long, npy_intp);

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
    char *use_norm          = args[11];
    char *use_fwd           = args[12];
    char *T_out             = args[13];

    npy_intp T_in_steps     = steps[0];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[13];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_octic;}
    else {FresT = &_fresnel_octic_norm;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    if (*(int *)use_fwd == 1){
        for (i=0; i<=*(long *)n_used; ++i){
            *(double *)(kd_vals+i*kd_steps) *= -1.0;
        }
    }

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
        P_4 = (4.375*P12-7.5)*P12+0.375;
        P_5 = ((7.875*P12-8.75)*P12+1.875)*P_1;
        P_6 = (((14.4375*P12-19.6875)*P12+6.5625)*P12)-0.3125;
        P_7 = (((P12*26.8125-43.3125)*P12+19.6875)*P12-2.1875)*P_1;

        /* Second set of polynomials, (P_{n} - P_1 * P_{n+1}) / (n+2) */
        b_0 = 0.5 - 0.5*P12;
        b_1 = (0.333333333333 - 0.333333333333*P_2)*P_1;
        b_2 = (P_2 - P_1*P_3)*0.25;
        b_3 = (P_3 - P_1*P_4)*0.20;
        b_4 = (P_4 - P_1*P_5)*0.16666666666666666;
        b_5 = (P_5 - P_1*P_6)*0.14285714285714285;
        b_6 = (P_6 - P_1*P_7)*0.125;

        /* Compute coefficients for Fresnel-Legendre Expansion. */
        A_0 = b_0 - Legendre_Coeff*P12;
        A_1 = b_1 - Legendre_Coeff*2.0*P_1*P_2;
        A_2 = b_2 - Legendre_Coeff*P_2*P_2;
        A_3 = b_3 - Legendre_Coeff*2.0*(P_1*P_4+P_2*P_3);
        A_4 = b_4 - Legendre_Coeff*(2.0*P_2*P_4+P_3*P_3);
        A_5 = b_5 - Legendre_Coeff*2.0*P_3*P_4;
        A_6 = b_6 - Legendre_Coeff*P_4*P_4;

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
        *((complex double *)T_out) = FresT(x_arr, T_in, w_func, rcpr_D,
                                           A_0, A_1, A_2, A_3, A_4, A_5, A_6,
                                           dx, rcpr_F,
                                           *(double *)kd_vals,
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
    free(x_arr);
    free(w_func);
}

static void Fresnel_Transform_Newton_Func(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data)
{
    long i, j, nw_pts, toler;
    double w_init, dx, two_dx, rcpr_F, EPS;

    toler = 5;
    EPS = 1.E-4;

    double (*fw)(double, double);
    complex double (*FresT)(double*, double*, char*, double*, double,
                            double, double, double, double, long, double,
                            double, long, npy_intp);

    char *T_in              = args[0];
    char *rho_km_vals       = args[1];
    char *F_km_vals         = args[2];
    char *phi_rad_vals      = args[3];
    char *kd_vals           = args[4];
    char *B_rad_vals        = args[5];
    char *D_km_vals         = args[6];
    char *w_km_vals         = args[7];
    char *start             = args[8];
    char *n_used            = args[9];
    char *wtype             = args[10];
    char *use_norm          = args[11];
    char *use_fwd           = args[12];
    char *T_out             = args[13];

    npy_intp T_in_steps     = steps[0];
    npy_intp rho_steps      = steps[1];
    npy_intp F_steps        = steps[2];
    npy_intp phi_steps      = steps[3];
    npy_intp kd_steps       = steps[4];
    npy_intp B_steps        = steps[5];
    npy_intp D_steps        = steps[6];
    npy_intp w_steps        = steps[7];
    npy_intp T_out_steps    = steps[13];

    if (*(int *)wtype == 0){fw = &__rect;}
    else if (*(int *)wtype == 1){fw = &__coss;}
    else if (*(int *)wtype == 2){fw = &Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 3){fw = &Kaiser_Bessel_Window_2_5;}
    else if (*(int *)wtype == 4){fw = &Kaiser_Bessel_Window_3_5;}
    else if (*(int *)wtype == 5){fw = &Modified_Kaiser_Bessel_Window_2_0;}
    else if (*(int *)wtype == 6){fw = &Modified_Kaiser_Bessel_Window_2_5;}
    else {fw = &Modified_Kaiser_Bessel_Window_3_5;}

    if (*(int *)use_norm == 0){FresT = &_fresnel_transform_newton;}
    else {FresT = &_fresnel_transform_newton_norm;}

    /* Compute first window width and window function. */
    phi_rad_vals    += *(long *)start * phi_steps;
    rho_km_vals     += *(long *)start * rho_steps;
    kd_vals         += *(long *)start * kd_steps;
    B_rad_vals      += *(long *)start * B_steps;
    D_km_vals       += *(long *)start * D_steps;
    F_km_vals       += *(long *)start * F_steps;
    w_km_vals       += *(long *)start * w_steps;
    T_in            += *(long *)start * T_in_steps;
    T_out           += *(long *)start * T_out_steps;

    if (*(int *)use_fwd == 1){
        for (i=0; i<=*(long *)n_used; ++i){
            *(double *)(kd_vals+i*kd_steps) *= -1.0;
        }
    }

    w_init  = *(double *)w_km_vals;
    dx      = *(double *)(rho_km_vals+rho_steps) - *(double *)(rho_km_vals);
    two_dx  = 2.0*dx;
    nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

    double* x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double* phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double* w_func  = (double *)malloc(sizeof(double) * nw_pts);

    for (j=0; j<nw_pts; ++j){
        x_arr[j]   =  *(double *)(rho_km_vals+rho_steps*(j-(nw_pts-1)/2));
        phi_arr[j] =  *(double *)(phi_rad_vals+phi_steps*(j-(nw_pts-1)/2));
        w_func[j]  = fw(x_arr[j] - *(double *)rho_km_vals, w_init);
    }

    for (i=0; i<=*(long *)n_used; ++i){
        rcpr_F          = 1.0 / *(double *)F_km_vals;

        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - *(double *)w_km_vals) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = *(double *)w_km_vals;
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;
            w_func  = (double *)realloc(w_func, sizeof(double)*nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            for (j=0; j<nw_pts; ++j){
                x_arr[j]  =  *(double *)(rho_km_vals+rho_steps*(j-(nw_pts-1)/2));
                phi_arr[j] =  *(double *)(phi_rad_vals+phi_steps*(j-(nw_pts-1)/2));
                w_func[j] = fw(x_arr[j] - *(double *)rho_km_vals, w_init);
            }
        }
        else {
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   =  *(double *)(rho_km_vals+rho_steps*(j-(nw_pts-1)/2));
                phi_arr[j] =  *(double *)(phi_rad_vals+phi_steps*(j-(nw_pts-1)/2));
                w_func[j]  = fw(x_arr[j] - *(double *)rho_km_vals, w_init);
            }
        }

        /*  Compute the fresnel tranform about the current point.   */
        *((complex double *)T_out) = FresT(x_arr, phi_arr, T_in, w_func,
                                           *(double *)kd_vals,
                                           *(double *)rho_km_vals,
                                           *(double *)B_rad_vals,
                                           *(double *)D_km_vals, EPS, toler,
                                           dx, rcpr_F, nw_pts, T_in_steps);

        /*  Increment pointers using pointer arithmetic, equivalent to        *
         *  changing var[n] to var[n+1].                                      */
        phi_rad_vals    += phi_steps;
        rho_km_vals     += rho_steps;
        kd_vals         += kd_steps;
        B_rad_vals      += B_steps;
        D_km_vals       += D_steps;
        F_km_vals       += F_steps;
        w_km_vals       += w_steps;
        T_in            += T_in_steps;
        T_out           += T_out_steps;
    }
    free(x_arr);
    free(w_func);
}

/*                            C-Python API Stuff                              */
static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

/* Define pointers to the C functions. */
PyUFuncGenericFunction f_quad_funcs[1] = {&Fresnel_Transform_Quadratic_Func};
PyUFuncGenericFunction f_cubic_funcs[1] = {&Fresnel_Transform_Cubic_Func};
PyUFuncGenericFunction f_quart_funcs[1] = {&Fresnel_Transform_Quartic_Func};
PyUFuncGenericFunction f_sxtic_funcs[1] = {&Fresnel_Transform_Sextic_Func};
PyUFuncGenericFunction f_octic_funcs[1] = {&Fresnel_Transform_Octic_Func};
PyUFuncGenericFunction f_newtn_funcs[1] = {&Fresnel_Transform_Newton_Func};

/* Input and return types for Quadratic Fresnel Transform */
static char quad_data_types[10] = {
    NPY_COMPLEX128,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_LONG,
    NPY_COMPLEX128
};

/* Input and return types for Quartic Fresnel Transform */
static char legendre_data_types[14] = {
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
    PyObject *fresnel_transform_cubic;
    PyObject *fresnel_transform_sextic;
    PyObject *fresnel_transform_octic;
    PyObject *fresnel_transform_newton;

    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_transform_quadratic = PyUFunc_FromFuncAndData(
        f_quad_funcs, PyuFunc_data, quad_data_types,
        1, 9, 1, PyUFunc_None, "fresnel_transform_quadratic",
        "fresnel_transform_quadratic_docstring", 0
    );

    fresnel_transform_cubic = PyUFunc_FromFuncAndData(
        f_cubic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_cubic",
        "fresnel_transform_cubic_docstring", 0
    );

    fresnel_transform_quartic = PyUFunc_FromFuncAndData(
        f_quart_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_quartic",
        "fresnel_transform_quartic_docstring", 0
    );

    fresnel_transform_sextic = PyUFunc_FromFuncAndData(
        f_sxtic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_sextic",
        "fresnel_transform_sextic_docstring", 0
    );

    fresnel_transform_octic = PyUFunc_FromFuncAndData(
        f_octic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_octic",
        "fresnel_transform_octic_docstring", 0
    );

    fresnel_transform_newton = PyUFunc_FromFuncAndData(
        f_newtn_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_newton",
        "fresnel_transform_newton_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform_quadratic", fresnel_transform_quadratic);
    PyDict_SetItemString(d, "fresnel_transform_cubic", fresnel_transform_cubic);
    PyDict_SetItemString(d, "fresnel_transform_quartic", fresnel_transform_quartic);
    PyDict_SetItemString(d, "fresnel_transform_sextic", fresnel_transform_sextic);
    PyDict_SetItemString(d, "fresnel_transform_octic", fresnel_transform_octic);
    PyDict_SetItemString(d, "fresnel_transform_newton", fresnel_transform_newton);
    Py_DECREF(fresnel_transform_quadratic);
    Py_DECREF(fresnel_transform_cubic);
    Py_DECREF(fresnel_transform_quartic);
    Py_DECREF(fresnel_transform_sextic);
    Py_DECREF(fresnel_transform_octic);
    Py_DECREF(fresnel_transform_newton);
    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_transform_quadratic;
    PyObject *fresnel_transform_quartic;
    PyObject *fresnel_transform_cubic;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _diffraction_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_transform_quadratic = PyUFunc_FromFuncAndData(
        f_quad_funcs, PyuFunc_data, quad_data_types,
        1, 9, 1, PyUFunc_None, "fresnel_transform_quadratic",
        "fresnel_transform_quadratic_docstring", 0
    );

    fresnel_transform_cubic = PyUFunc_FromFuncAndData(
        f_cubic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_cubic",
        "fresnel_transform_cubic_docstring", 0
    );

    fresnel_transform_quartic = PyUFunc_FromFuncAndData(
        f_quart_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_quartic",
        "fresnel_transform_quartic_docstring", 0
    );

    fresnel_transform_sextic = PyUFunc_FromFuncAndData(
        f_sxtic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_sextic",
        "fresnel_transform_sextic_docstring", 0
    );

    fresnel_transform_octic = PyUFunc_FromFuncAndData(
        f_octic_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_octic",
        "fresnel_transform_octic_docstring", 0
    );

    fresnel_transform_octic = PyUFunc_FromFuncAndData(
        f_newtn_funcs, PyuFunc_data, legendre_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_newton",
        "fresnel_transform_newton_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform_quadratic", fresnel_transform_quadratic);
    PyDict_SetItemString(d, "fresnel_transform_cubic", fresnel_transform_cubic);
    PyDict_SetItemString(d, "fresnel_transform_quartic", fresnel_transform_quartic);
    PyDict_SetItemString(d, "fresnel_transform_sextic", fresnel_transform_sextic);
    PyDict_SetItemString(d, "fresnel_transform_octic", fresnel_transform_octic);
    PyDict_SetItemString(d, "fresnel_transform_newton", fresnel_transform_newton);
    Py_DECREF(fresnel_transform_quadratic);
    Py_DECREF(fresnel_transform_cubic);
    Py_DECREF(fresnel_transform_quartic);
    Py_DECREF(fresnel_transform_sextic);
    Py_DECREF(fresnel_transform_octic);
    Py_DECREF(fresnel_transform_newton);
    return m;
}
#endif
