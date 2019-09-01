/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H
#define RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H

/*  Various trig functions, complex variables, and more found here.           */
#include <math.h>
#include <complex.h>

/*  Coefficients and constants defined here.                                  */
#include "__math_constants.h"

/*  Functions for computing the Fresnel Kernel and Newton's Method.           */
#include "__fresnel_kernel.h"

extern complex double Fresnel_Transform_Double(
    double *x_arr, char *T_in, double *w_func, double F, double dx,
    long n_pts, long T_in_steps
);

extern complex double Fresnel_Transform_Norm_Double(
    double *x_arr, char *T_in, double *w_func, double F,
    double dx, long n_pts, long T_in_steps
);

extern complex double Fresnel_Legendre_Double(
    double *x_arr, char *T_in, double *w_func, double D, double *coeffs,
    double dx, double F, double kd, long n_pts, int order, long T_in_steps
);

extern complex double Fresnel_Legendre_Norm_Double(
    double *x_arr, char *T_in, double *w_func, double D, double *coeffs,
    double dx, double F, double kd, long n_pts, int order, long T_in_steps
);

extern complex double Fresnel_Transform_Newton_Double(
    double *x_arr, double *phi_arr, char *T_in, double *w_func, double kD,
    double r, double B, double D, double EPS, long toler, double dx, double F,
    long n_pts, long T_in_steps
);

extern complex double Fresnel_Transform_Newton_Norm_Double(
    double *x_arr, double *phi_arr, char *T_in, double *w_func, double kD,
    double r, double B, double D, double EPS, long toler, double dx, double F,
    long n_pts, long T_in_steps
);


#endif