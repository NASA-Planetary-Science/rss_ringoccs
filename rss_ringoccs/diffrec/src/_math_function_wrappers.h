
/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = BesselJ0_Float(in[i]);
    }
}

static void double_besselJ0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = BesselJ0_Double(in[i]);
    }
}

static void long_double_besselJ0(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data){
    npy_intp i;
    npy_intp n = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = BesselJ0_Long_Double(in[i]);
    }
}