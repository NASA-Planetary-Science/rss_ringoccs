#ifndef RSS_RINGOCCS_SINC_WRAPPERS_H
#define RSS_RINGOCCS_SINC_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_sinc(char **args, npy_intp *dimensions,
                       npy_intp *steps, void *data)
{
    npy_intp dim = dimensions[0];

    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    Get_Float_Array(in, out, dim, Sinc_Float);
}

static void double_sinc(char **args, npy_intp *dimensions,
                        npy_intp *steps, void *data)
{
    npy_intp dim = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    Get_Double_Array(in, out, dim, Sinc_Double);
}

static void long_double_sinc(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data)
{
    npy_intp dim = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    Get_Long_Double_Array(in, out, dim, Sinc_Long_Double);
}

#endif