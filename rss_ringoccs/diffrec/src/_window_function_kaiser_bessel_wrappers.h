/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_WINDOW_FUNCTION_KAISER_BESSEL_WRAPPERS_H
#define RSS_RINGOCCS_WINDOW_FUNCTION_KAISER_BESSEL_WRAPPERS_H

#include "__window_functions.h"

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

static void float_kb20(char **args, npy_intp *dimensions,
                       npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *x   =  (float *)args[0];
    float  W   = *(float *)args[1];
    float *out =  (float *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_0_Float(x[i], W);
    }
}

static void double_kb20(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_0_Double(x[i], W);
    }
}

static void long_double_kb20(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *x   =  (long double *)args[0];
    long double  W   = *(long double *)args[1];
    long double *out =  (long double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_0_Long_Double(x[i], W);
    }
}

static void float_kb25(char **args, npy_intp *dimensions,
                       npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *x   =  (float *)args[0];
    float  W   = *(float *)args[1];
    float *out =  (float *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_5_Float(x[i], W);
    }
}

static void double_kb25(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_5_Double(x[i], W);
    }
}

static void long_double_kb25(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *x   =  (long double *)args[0];
    long double  W   = *(long double *)args[1];
    long double *out =  (long double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_2_5_Long_Double(x[i], W);
    }
}

static void float_kb35(char **args, npy_intp *dimensions,
                       npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *x   =  (float *)args[0];
    float  W   = *(float *)args[1];
    float *out =  (float *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_3_5_Float(x[i], W);
    }
}

static void double_kb35(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x   =  (double *)args[0];
    double  W   = *(double *)args[1];
    double *out =  (double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_3_5_Double(x[i], W);
    }
}

static void long_double_kb35(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *x   =  (long double *)args[0];
    long double  W   = *(long double *)args[1];
    long double *out =  (long double *)args[2];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_3_5_Long_Double(x[i], W);
    }
}

static void float_kbal(char **args, npy_intp *dimensions,
                       npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *x    =  (float *)args[0];
    float  W    = *(float *)args[1];
    float alpha = *(float *)args[2];
    float *out  =  (float *)args[3];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_Al_Float(x[i], W, alpha);
    }
}

static void double_kbal(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *x    =  (double *)args[0];
    double  W    = *(double *)args[1];
    double alpha = *(double *)args[2];
    double *out  =  (double *)args[3];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_Al_Double(x[i], W, alpha);
    }
}

static void long_double_kbal(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *x    =  (long double *)args[0];
    long double  W    = *(long double *)args[1];
    long double alpha = *(long double *)args[2];
    long double *out  =  (long double *)args[3];

    for (i = 0; i < n; i++) {
        out[i] = Kaiser_Bessel_Al_Long_Double(x[i], W, alpha);
    }
}

#endif