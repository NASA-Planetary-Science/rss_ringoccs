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

static void __rect(char *wfunc, double w_width, double dx, npy_intp w_steps)
{
        int i, nw_pts;
        double x;

        /* Window functions have an odd number of points. */
        nw_pts = (2 * (int)(w_width / (2.0 * dx))) + 1;

        for (i=0; i<nw_pts; i++){
            *((double *)(wfunc+i*w_steps) = 1.0;
        }
}

static void __coss(char *wfunc, double w_width, double dx,
                   npy_intp w_steps, char *nw_pts)
{
        int i;
        double x;

        /* Window functions have an odd number of points. */
        *(int *)nw_pts = (2 * (int)(w_width / (2.0 * dx))) + 1;

        for (i=0; i<nw_pts; i++){
            x = i - (nw_pts - 1) / 2.0);
            x *= ONE_PI * dx / w_width;
            x = cos(x);
            *((double *)(wfunc+i*w_steps) = x*x;
        }
}

static PyMethodDef _special_functions_methods[] = {{NULL, NULL, 0, NULL}};

static void Fresnel_Transform_Func(char **args, npy_intp *dimensions,
                                      npy_intp* steps, void* data)
{
    npy_intp i, n_total;
    npy_intp n = dimensions[0];

    char *T_in          = args[0];
    char *rho_km_vals   = args[1];
    char *F_km_vals     = args[2];
    char *w_km_vals     = args[3];
    char *w_func        = args[4];
    char *x_arr         = args[5];
    char *start         = args[6];
    char *n_used        = args[7];
    char *dx_km         = args[8];
    char *T_out         = args[9];

    npy_intp T_in_steps     = steps[0];
    npy_intp rho_steps      = steps[1];
    npy_intp F_steps        = steps[2];
    npy_intp w_steps        = steps[3];
    npy_intp T_out_steps    = steps[7];

    char *n_pts;
    double kD_vals, w_init, w_width, r, F;
    static void (*fw)(char *, double, double, npy_intp);

    n_total = n_used*rho_steps;

    if (wtype == "rect"){
        fw = &__rect;
    }
    else if (strcomp(wtype, "coss", 4) == 0){
        fw = &__coss;
    }
    else if (strcomp(wtype, "kb20", 4) == 0){
        fw = &__coss;
    }
    else if (strcomp(wtype, "kb25", 4) == 0){
        fw = &__coss;
    }
    else if (strcomp(wtype, "kb35", 4) == 0){
        fw = &__coss;
    }
    else if (strcomp(wtype, "kbmd20", 4) == 0){
        fw = &__coss;
    }
    else {
        fw = &__coss;
    }

    /* Compute first window width and window function. */
    w_km_vals += start*w_steps;
    rho_km_vals += start*rho_steps;
    F_km_vals += start*F_steps;
    w_init = *(double *)w_km_vals;
    fw(*w_func, w_init, *(double *)dx_km, w_steps, *n_pts);

    for (i=0; i<=n_total; ++i){
        /* Rho, Window width, Frensel scale for current point. */
        r       = *(double *)rho_km_vals;
        w_width = *(double *)w_km_vals;
        F       = *(double *)F_km_vals;

        if (fabs(w_init - w_width) >= 2.0 * *(double *)dx_km) {
            /* Reset w_init and recompute window function. */
            w_init = w_width;
            fw(*w_func, w_init, *(double *)dx_km, w_steps, *n_pts);

            /* Compute psi for with stationary phase value */
            x = r-r0
            x2 = HALF_PI * x * x
        }
        else {
            crange += 1   
            psi_vals = x2 / F2[center]
        }

        /* Compute kernel function for Fresnel inverse. */
        ker = w_func*np.exp(-1j*psi_vals)

        /* Compute 'approximate' Fresnel Inversion for current point. */
        T_out[center] = np.sum(ker*T)*self.dx_km*(0.5+0.5j)/F

        rho_km_vals += rho_steps;
        w_km_vals   += w_steps;
        F_km_vals   += F_steps;
    }