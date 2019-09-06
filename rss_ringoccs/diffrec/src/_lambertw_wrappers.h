/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_LAMBERTW_WRAPPERS_H
#define RSS_RINGOCCS_LAMBERTW_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_lambertw(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, LambertW_Float);
}

static void double_lambertw(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements, LambertW_Double);
}

static void long_double_lambertw(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements, LambertW_Long_Double);
}

#endif