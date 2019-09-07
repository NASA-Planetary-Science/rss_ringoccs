/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RESOLUTION_INVERSE_FUNCTION_WRAPPERS_H
#define RSS_RINGOCCS_RESOLUTION_INVERSE_FUNCTION_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_resolution_inverse(char **args, npy_intp *dimensions,
                                     npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, Resolution_Inverse_Float);
}

static void double_resolution_inverse(char **args, npy_intp *dimensions,
                                      npy_intp *steps, void *data)
{
    npy_intp dim = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    Get_Double_Array(in, out, dim, Resolution_Inverse_Double);
}

static void long_double_resolution_inverse(char **args, npy_intp *dimensions,
                                           npy_intp *steps, void *data)
{
    npy_intp dim = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    Get_Long_Double_Array(in, out, dim, Resolution_Inverse_Long_Double);
}

#endif