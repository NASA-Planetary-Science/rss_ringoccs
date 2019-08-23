/******************************************************************************
 *                          Diffraction Functions                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains functions used for computing the Fresnel Inverse   *
 *      Transform on a set of diffraction limited data. There are several     *
 *      Methods of performing this:                                           *
 *          Fresnel Quadratic Approximations:                                 *
 *              Classic quadratic approximation from Fourier Optics.          *
 *          Legendre Cubic Expansion:                                         *
 *              Approximation of the Fresnel kernel by a quartic polynomial   *
 *              using the Legendre Polynomials.                               *
 *          Legendre Quartic Expansion:                                       *
 *              Approximation of the Fresnel kernel by a quartic polynomial   *
 *              using the Legendre Polynomials.                               *
 *          Legendre Sextic Expansion:                                        *
 *              Approximation of the Fresnel kernel by a sextic polynomial    *
 *              using the Legendre Polynomials.                               *
 *          Legendre Octic Expansion:                                         *
 *              Approximation of the Fresnel kernel by a octic polynomial     *
 *              using the Legendre Polynomials.                               *
 ******************************************************************************
 *  Variables:                                                                *
 *     A_0         (Double):                                                  *
 *         The coefficient of the x^2 term in the expansion for psi.          *
 *     A_1         (Double):                                                  *
 *         The coefficient of the x^3 term in the expansion for psi.          *
 *     A_2         (Double):                                                  *
 *         The coefficient of the x^4 term in the expansion for psi.          *
 *     A_3         (Double):                                                  *
 *         The coefficient of the x^5 term in the expansion for psi.          *
 *     A_4         (Double):                                                  *
 *         The coefficient of the x^6 term in the expansion for psi.          *
 *     A_5         (Double):                                                  *
 *         The coefficient of the x^7 term in the expansion for psi.          *
 *     A_6         (Double):                                                  *
 *         The coefficient of the x^8 term in the expansion for psi.          *
 *     dx          (Double):                                                  *
 *         Spacing between two consecutive ring intercept points, km.         *
 *     kd          (Double):                                                  *
 *         The wavenumber, k, weighted by the spacecraft-ring distance, D.    *
 *     n_pts       (Long):                                                    *
 *         Half the number of points in the window, rounded down.             *
 *     rcpr_D      (Double):                                                  *
 *         1/D, where D is the distant from the spacecraft to the ring        *
 *         intercept point, in kilometers.                                    *
 *     rcpr_F      (Double):                                                  *
 *         The reciprocal of the Fresnel scale in kilometers.                 *
 *     T_in        (Pointer to Char):                                         *
 *         The raw diffraction-limited data that is to be corrected.          *
 *     T_in_steps  (npy_intp (Equivalent to Long)):                           *
 *         The number of steps in memory to get from the nth data point       *
 *         to the (n+1)th data point in the T_in variable (See above).        *
 *     T_out       (Complex Double):                                          *
 *         The diffraction corrected profile.                                 *
 *     w_func      (Pointer to Double):                                       *
 *         Pre-computed window function. Should have the same number of       *
 *         points as the x_arr pointer.                                       *
 *     x_arr       (Pointer to Double):                                       *
 *         The ring radii within the given window.                            *
 ******************************************************************************
 *  The Inverse Fresnel Transform:                                            *
 *                                                                            *
 *                W/2                                                         *
 *                 -                                                          *
 *                | |                                                         *
 *     T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0               *
 *              | |                                                           *
 *               -                                                            *
 *              -W/2                                                          *
 *                                                                            *
 *  Where T_hat is the diffracted data, w is the window function, r is        *
 *  the ring intercept point, and r_0 is a dummy variable of integration.     *
 *  psi is the Fresnel Kernel, and exp is simply the exponential function.    *
 ******************************************************************************
 *  The Normalization Scheme:                                                 *
 *      As the resolution get's too high, say 10 km or greater, then window   *
 *      width quickly shrinks to zero. Thus the integral will be approximately*
 *      zero. To account for this, the option to normalize the integral by the*
 *      window width is offered. The normalization is defined as follows:     *
 *                                                                            *
 *                    |     _ +infinity           |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -  -infinity            |                           *
 *          Norm =  __________________________________                        *
 *                  |    -  +W/2                    |                         *
 *                  |   | |                         |                         *
 *                  |   |    w(x) exp(-i psi(x)) dx |                         *
 *                  | | |                           |                         *
 *                  |  -   -W/2                     |                         *
 *                                                                            *
 *      This has the effect of making the free-space regions, or regions which*
 *      were not affected by diffraction, evaluate to approximately one,      *
 *      regardless of what resolution was chosen.                             *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  get_arr:                                                                  *
 *      Void function that takes in a pointer to a double array and creates   *
 *      an array of values for the ring radius within half a window width of  *
 *      the ring intercept point. Do to symmetry, only the values to the left *
 *      of the ring intercept point are computed.                             *
 ******************************************************************************
 *  _fresnel_transform:                                                       *
 *      Computes the Fresnel Inverse Transform using the classic Fresnel      *
 *      quadratic approximation. No normalization is applied.                 *
 ******************************************************************************
 *  _fresnel_transform_norm:                                                  *
 *      Same as _fresnel_transform, but the normalization is applied.         *
 ******************************************************************************
 *  _fresnel_cubic,                                                           *
 *  _fresnel_quartic,                                                         *
 *  _fresnel_sextic,                                                          *
 *  _fresnel_octic:                                                           *
 *      Computes the Fresnel Inverse Transform using Legendre Polynomials to  *
 *      approximate the Fresnel kernel to various powers (3, 4, 6, or 8).     *
 ******************************************************************************
 *  _fresnel_cubic_norm,                                                      *
 *  _fresnel_quartic_norm,                                                    *
 *  _fresnel_sextic_norm,                                                     *
 *  _fresnel_octic_norm:                                                      *
 *      Same as previous functions, but with the normalization scheme.        *
 ******************************************************************************
 *                             A FRIENDY WARNING                              *
 ******************************************************************************
 *  This code uses complex numbers throughout, and is compatible with the C99 *
 *  standard. To use this code, make sure your compiler supports C99 or more  *
 *  recent standards of the C Programming Language.                           *
 *****************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  Various trig functions, complex variables, and more found here.          */
#include <math.h>
#include <complex.h>

/*  Various header files required for the C-Python API to work.              */
#include "../../../include/Python.h"
#include "../../../include/ndarraytypes.h"
#include "../../../include/ufuncobject.h"

/*  Various coefficients and constants defined here.                         */
#include "__math_constants.h"

/*  Window functions and Fresnel Transforms defined here.                    */
#include "__window_functions.h"
#include "__diffraction_functions.h"

/*  Functions for computing the Fresnel Kernel and Newton's Method.          */
#include "__fresnel_kernel.h"

static void Fresnel_Transform_Quadratic_Func(char **args, npy_intp *dimensions,
                                             npy_intp* steps, void* data)
{
    /*  i and j used for indexing, nw_pts is number of points in window.     */
    long i, j, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).  */
    double w_init, two_dx;

    /*  Declare function pointers. fw is window function, FresT is Fresnel   *
     *  transform. Options for FresT are with or without normalization.      */
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double,
                            double, long, npy_intp);

    /*  Grab all data passed in from Python.                                 */
    char *T_in            = args[0];
    double dx             = *(double *)args[1];
    double *F_km_vals     =  (double *)args[2];
    double *w_km_vals     =  (double *)args[3];
    long start            = *(long *)args[4];
    long n_used           = *(long *)args[5];
    int wtype             = *(int *)args[6];
    int use_norm          = *(int *)args[7];
    int use_fwd           = *(int *)args[8];
    complex double *T_out = (complex double *)args[9];

    /*  Number of steps in memory from the ith point to the (i+1)th point.   */
    npy_intp T_in_steps  = steps[0];

    /*  Cast the selected window type to the fw pointer.                     */
    if      (wtype == 0){fw = &Rect_Window_Double;}
    else if (wtype == 1){fw = &Coss_Window_Double;}
    else if (wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else    {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    /*  Cast FresT to the appropriate function.                              */
    if (use_norm == 0){FresT = &_fresnel_transform;}
    else {FresT = &_fresnel_transform_norm;}

    /* Move the pointers to the correct starting point. Compute window width.*/
    T_in   += start*T_in_steps;
    w_init  = w_km_vals[start];
    center  = start;

    /*  Compute some extra necessary variables.                              */
    two_dx = 2.0*dx;
    nw_pts = (long)(w_init / (two_dx));

    /*  Reserve some memory for two arrays, the ring radius and the window   *
     *  function. This will need to be reallocated later if the window width *
     *  changes by more than two_dx.                                         */
    double* x_arr  = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);

    /*  Pass the x_arr array (ring radius) to the void function get_arr.     *
     *  This alters the x_arr pointer so that it's values range from -W/2 to *
     *  zero, where W is the window width.                                   */
    get_arr(x_arr, dx, nw_pts);

    /* Compute Window Functions, and compute pi/2 * x^2                      */
    for(j=0; j<nw_pts; ++j){
        w_func[j] = fw(x_arr[j], w_init);
        x_arr[j] *= PI_BY_TWO*x_arr[j];

        /*  If forward transform is selected, negate x_arr.                  */
        x_arr[j] *= (use_fwd == 0) - (use_fwd == 1);
    }

    /*  Compute the Fresnel transform across the input data.                 */
    for (i=0; i<=n_used; ++i){
        if (fabs(w_init - w_km_vals[center]) >= two_dx) {
            /* Reset w_init and recompute window function.                   */
            w_init = w_km_vals[center];
            nw_pts = (long)(w_init / two_dx);

            /*  Reallocate memory, since the sizes of the arrays changed.    */
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.        */
            get_arr(x_arr, dx, nw_pts);

            /* Compute Window Functions, and compute pi/2 * x^2              */
            for(j=0; j<nw_pts; ++j){
                w_func[j] = fw(x_arr[j], w_init);
                x_arr[j] *= PI_BY_TWO*x_arr[j];
                x_arr[j] *= (use_fwd == 0) - (use_fwd == 1);
            }
        }

        /*  Compute the Fresnel Transform about the current point.           */
        T_out[center] = FresT(x_arr, T_in, w_func, F_km_vals[center], dx,
                              nw_pts, T_in_steps);

        /*  Move the pointers to the next point.                             */
        center += 1;
        T_in   += T_in_steps;
    }
    free(x_arr);
    free(w_func);
}

static void Fresnel_Legendre_Transform_Func(char **args, npy_intp *dimensions,
                                            npy_intp* steps, void* data)
{
    long i, j, nw_pts, center;
    int k, order_by_2;
    double w_init, two_dx, cosb, sinp, cosp;
    double Legendre_Coeff, rcpr_D, rcpr_F;
    double (*fw)(double, double);
    complex double (*FresT)(double*, char*, double*, double, double *,
                            double, double, double, long, int, npy_intp);

    char *T_in            = args[0];
    double dx             = *(double *)args[1];
    double *F_km_vals     =  (double *)args[2];
    double *phi_rad_vals  =  (double *)args[3];
    double *kd_vals       =  (double *)args[4];
    double *B_rad_vals    =  (double *)args[5];
    double *D_km_vals     =  (double *)args[6];
    double *w_km_vals     =  (double *)args[7];
    long start            = *(long *)args[8];
    long n_used           = *(long *)args[9];
    int wtype             = *(int *)args[10];
    int use_norm          = *(int *)args[11];
    int use_fwd           = *(int *)args[12];
    int order             = *(int *)args[13];
    complex double *T_out =  (complex double *)args[14];

    npy_intp T_in_steps     = steps[0];

    /*  Cast the selected window type to the fw pointer.                      */
    if      (wtype == 0){fw = &Rect_Window_Double;}
    else if (wtype == 1){fw = &Coss_Window_Double;}
    else if (wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else                {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    if (use_norm){FresT = &_fresnel_legendre_norm;}
    else         {FresT = &_fresnel_legendre;}

    /* Compute first window width and window function. */
    T_in  += start * T_in_steps;
    center = start;

    if (use_fwd){
        for (i=0; i<= n_used; ++i){
            kd_vals[i] *= -1.0;
        }
    }

    w_init = w_km_vals[center];
    two_dx = 2.0*dx;
    nw_pts = (long)(w_init / (two_dx));

    double* x_arr      = (double *)malloc(sizeof(double) * nw_pts);
    double* w_func     = (double *)malloc(sizeof(double) * nw_pts);
    double* legendre_p = (double *)malloc(sizeof(double) * (order+1));
    double* fresnel_p  = (double *)malloc(sizeof(double) * order);
    double* coeffs_p   = (double *)malloc(sizeof(double) * order);
    double l_coeffs[order];
    order_by_2 = (order+1)/2;

    for (j=0; j<order; ++j){
        l_coeffs[j] = 1.0/(j+2.0);
    }

    reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    legendre_p[0] = 1.0;
    for (i = 0; i <= n_used; ++i){
        rcpr_F          = 1.0 / F_km_vals[center];
        rcpr_D          = 1.0 / D_km_vals[center];
        cosb            = cos(B_rad_vals[center]);
        cosp            = cos(phi_rad_vals[center]);
        sinp            = sin(phi_rad_vals[center]);
        Legendre_Coeff  = cosb*sinp;
        Legendre_Coeff *= Legendre_Coeff;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        legendre_p[1] = cosb*cosp;
        fresnel_p[0]  = 0.5-0.5*legendre_p[1]*legendre_p[1];

        /* Compute Legendre Polynomials,                                      */
        for (j=1; j < order; ++j){
            legendre_p[j+1] = ((2.0*j+1.0)*legendre_p[1]*legendre_p[j] - j*legendre_p[j-1])*l_coeffs[j-1];
            fresnel_p[j] = (legendre_p[j]-legendre_p[1]*legendre_p[j+1]) * l_coeffs[j];
        }

        /*  Compute the coefficients using Cauchy Products. First compute     *
         *  the bottom triangle of the square in the product.                 */
        for (j=1; j<=order_by_2; ++j){
            coeffs_p[j-1] = 0.0;
            for (k=0; k<j; k++){
                coeffs_p[j-1] += legendre_p[k+1]*legendre_p[j-k];
            }
            coeffs_p[j-1] = fresnel_p[j-1] - Legendre_Coeff*coeffs_p[j-1];
        }

        /*  Compute along the upper triangle of the square.                   */
        for (j=order_by_2+1; j<order; ++j){
            coeffs_p[j-1] = 0.0;
            for (k=j-order_by_2; k<order_by_2; k++){
                coeffs_p[j-1] += legendre_p[k+1]*legendre_p[j-k];
            }
            coeffs_p[j-1] = fresnel_p[j-1] - Legendre_Coeff*coeffs_p[j-1];
        }

        /* Compute the last coefficient.                                      */
        coeffs_p[order-1] = legendre_p[order_by_2]*legendre_p[order_by_2];
        coeffs_p[order-1] = fresnel_p[order-1]-Legendre_Coeff*coeffs_p[order-1];


        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - w_km_vals[center]) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init = w_km_vals[center];
            nw_pts = (long)(w_init / two_dx);
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);
        }

        /*  Compute the fresnel tranform about the current point.   */
        T_out[center] = FresT(x_arr, T_in, w_func, rcpr_D, coeffs_p, dx, rcpr_F,
                              kd_vals[center], nw_pts, order, T_in_steps);

        /*  Increment pointers using pointer arithmetic.                     */
        T_in   += T_in_steps;
        center += 1;
    }

    free(x_arr);
    free(w_func);
    free(legendre_p);
    free(fresnel_p);
    free(coeffs_p);
}

static void Fresnel_Transform_Newton_Func(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data)
{
    long i, j, nw_pts, toler, center;
    double w_init, dx, two_dx, rcpr_F, EPS;

    toler = 5;
    EPS = 1.E-4;

    double (*fw)(double, double);
    complex double (*FresT)(double*, double*, char*, double*, double,
                            double, double, double, double, long, double,
                            double, long, npy_intp);

    char *T_in            = args[0];
    double *rho_km_vals   = (double *)args[1];
    double *F_km_vals     = (double *)args[2];
    double *phi_rad_vals  = (double *)args[3];
    double *kd_vals       = (double *)args[4];
    double *B_rad_vals    = (double *)args[5];
    double *D_km_vals     = (double *)args[6];
    double *w_km_vals     = (double *)args[7];
    long start            = *(long *)args[8];
    long n_used           = *(long *)args[9];
    int wtype             = *(int *)args[10];
    int use_norm          = *(int *)args[11];
    int use_fwd           = *(int *)args[12];
    complex double *T_out = (complex double *)args[13];

    npy_intp T_in_steps     = steps[0];

    /*  Cast the selected window type to the fw pointer.                      */
    if      (wtype == 0){fw = &Rect_Window_Double;}
    else if (wtype == 1){fw = &Coss_Window_Double;}
    else if (wtype == 2){fw = &Kaiser_Bessel_2_0_Double;}
    else if (wtype == 3){fw = &Kaiser_Bessel_2_5_Double;}
    else if (wtype == 4){fw = &Kaiser_Bessel_3_5_Double;}
    else if (wtype == 5){fw = &Modified_Kaiser_Bessel_2_0_Double;}
    else if (wtype == 6){fw = &Modified_Kaiser_Bessel_2_5_Double;}
    else                {fw = &Modified_Kaiser_Bessel_3_5_Double;}

    if (use_norm){FresT = &_fresnel_transform_newton_norm;}
    else {FresT = &_fresnel_transform_newton;}

    /* Compute first window width and window function. */
    center = start;
    T_in  += start * T_in_steps;

    if (use_fwd){
        for (i=0; i<=n_used; ++i){
            kd_vals[center+i] *= -1.0;
        }
    }

    w_init  = w_km_vals[center];
    dx      = rho_km_vals[center+1] - rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    double *x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double *phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double *w_func  = (double *)malloc(sizeof(double) * nw_pts);

    for (j=0; j<nw_pts; ++j){
        x_arr[j]   = rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - rho_km_vals[center], w_init);
    }

    for (i=0; i<=n_used; ++i){
        rcpr_F = 1.0 / F_km_vals[center];

        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - w_km_vals[center]) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = phi_rad_vals[center+j-(nw_pts-1)/2];
                w_func[j]  = fw(x_arr[j] - rho_km_vals[center], w_init);
            }
        }
        else {
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = phi_rad_vals[center+j-(nw_pts-1)/2];
            }
        }

        /*  Compute the fresnel tranform about the current point.   */
        T_out[center] = FresT(x_arr, phi_arr, T_in, w_func, kd_vals[center],
                              rho_km_vals[center], B_rad_vals[center],
                              D_km_vals[center], EPS, toler, dx, rcpr_F,
                              nw_pts, T_in_steps);

        /*  Increment pointers using pointer arithmetic.                      */
        T_in   += T_in_steps;
        center += 1;
    }
    free(x_arr);
    free(phi_arr);
    free(w_func);
}

/*                            C-Python API Stuff                              */
static PyMethodDef _diffraction_functions_methods[] = {{NULL, NULL, 0, NULL}};

/* Define pointers to the C functions. */
PyUFuncGenericFunction f_quad_funcs[1]     = {&Fresnel_Transform_Quadratic_Func};
PyUFuncGenericFunction f_legendre_funcs[1] = {&Fresnel_Legendre_Transform_Func};
PyUFuncGenericFunction f_newtn_funcs[1]    = {&Fresnel_Transform_Newton_Func};

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
static char legendre_data_types[15] = {
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
    NPY_LONG,
    NPY_COMPLEX128
};

static char newton_data_types[14] = {
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
    PyObject *fresnel_legendre_transform;
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

    fresnel_legendre_transform = PyUFunc_FromFuncAndData(
        f_legendre_funcs, PyuFunc_data, legendre_data_types,
        1, 14, 1, PyUFunc_None, "fresnel_transform_cubic",
        "fresnel_transform_cubic_docstring", 0
    );


    fresnel_transform_newton = PyUFunc_FromFuncAndData(
        f_newtn_funcs, PyuFunc_data, newton_data_types,
        1, 13, 1, PyUFunc_None, "fresnel_transform_newton",
        "fresnel_transform_newton_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_transform_quadratic", fresnel_transform_quadratic);
    PyDict_SetItemString(d, "fresnel_legendre_transform", fresnel_legendre_transform);
    PyDict_SetItemString(d, "fresnel_transform_newton", fresnel_transform_newton);
    Py_DECREF(fresnel_transform_quadratic);
    Py_DECREF(fresnel_legendre_transform);
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
