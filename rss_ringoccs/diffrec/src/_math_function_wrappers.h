
/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>
#include "__math_functions.h"

static void float_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, BesselJ0_Float);
}

static void double_besselJ0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{

    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements, BesselJ0_Double);
}

static void long_double_besselJ0(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements, BesselJ0_Long_Double);
}

static void float_besselI0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, BesselI0_Float);
}

static void double_besselI0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements, BesselI0_Double);
}

static void long_double_besselI0(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements, BesselI0_Long_Double);
}