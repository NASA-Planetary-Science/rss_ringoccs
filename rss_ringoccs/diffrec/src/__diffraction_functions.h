/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H
#define RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H

/*  Various trig functions, complex variables, and more found here.           */
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "__math_functions.h"
#include "__window_functions.h"
#include "__fresnel_kernel.h"

typedef struct _dlpdataobj {
    complex double *T_in;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_rad_vals;
    double *kd_vals;
    double *B_rad_vals;
    double *D_km_vals;
    double *w_km_vals;
    long start;
    long n_used;
    int wtype;
    int use_norm;
    int use_fwd;
    int psitype;
    int order;
    complex double *T_out;
} DLPObj;

/*  Coefficients and constants defined here.                                  */
#include "__math_constants.h"

/*  Functions for computing the Fresnel Kernel and Newton's Method.           */
#include "__fresnel_kernel.h"

extern complex double Fresnel_Transform_Double(
    double *x_arr, complex double *T_in, double *w_func, double F, double dx,
    long n_pts, long center
);

extern complex double Fresnel_Transform_Norm_Double(
    double *x_arr, complex double *T_in, double *w_func, double F,
    double dx, long n_pts, long center
);

extern complex double Fresnel_Legendre_Double(
    double *x_arr, complex double *T_in, double *w_func, double D,
    double *coeffs, double dx, double F, double kd, long n_pts,
    int order, long center
);

extern complex double Fresnel_Legendre_Norm_Double(
    double *x_arr, complex double *T_in, double *w_func, double D,
    double *coeffs, double dx, double F, double kd, long n_pts,
    int order, long center
);

extern complex double Fresnel_Transform_Newton_Double(
    double *x_arr, double *phi_arr, complex double *T_in, double *w_func,
    double kD, double r, double B, double D, double EPS, long toler, double dx,
    double F, long n_pts, long center
);

extern complex double Fresnel_Transform_Newton_Norm_Double(
    double *x_arr, double *phi_arr, complex double *T_in, double *w_func,
    double kD, double r, double B, double D, double EPS, long toler, double dx,
    double F, long n_pts, long center
);

extern void DiffractionCorrectionFresnel(DLPObj dlp);

extern void DiffractionCorrectionLegendre(DLPObj dlp);

extern void DiffractionCorrectionNewton(DLPObj dlp);

#endif