/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RESOLUTION_INVERSE_FUNCTION_WRAPPERS_H
#define RSS_RINGOCCS_RESOLUTION_INVERSE_FUNCTION_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_resolution_inverse(char **args, npy_intp *dimensions,
                                     npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Resolution_Inverse_Float(in[i]);
    }
}

static void double_resolution_inverse(char **args, npy_intp *dimensions,
                                      npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Resolution_Inverse_Double(in[i]);
    }
}

static void long_double_resolution_inverse(char **args, npy_intp *dimensions,
                                           npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Resolution_Inverse_Long_Double(in[i]);
    }
}

#endif