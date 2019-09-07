#ifndef RSS_RINGOCCS_FRESNEL_INTEGRALS_WRAPPERS_H
#define RSS_RINGOCCS_FRESNEL_INTEGRALS_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

#include "__fresnel_integrals.h"

/*--------------------------Fresnel Sine Function-----------------------------*/

static void float_fresnelsin(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, Fresnel_Sine_Taylor_to_Asymptotic_Float);
}

static void double_fresnelsin(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements,
                     Fresnel_Sine_Taylor_to_Asymptotic_Double);
}

static void long_double_fresnelsin(char **args, npy_intp *dimensions,
                                   npy_intp* steps, void* data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];


    Get_Long_Double_Array(x, y, n_elements,
                          Fresnel_Sine_Taylor_to_Asymptotic_Long_Double);
}

/*-------------------------Fresnel Cosine Function----------------------------*/

static void float_fresnelcos(char **args, npy_intp *dimensions,
                             npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements,
                    Fresnel_Cosine_Taylor_to_Asymptotic_Float);
}

static void double_fresnelcos(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements,
                     Fresnel_Cosine_Taylor_to_Asymptotic_Double);
}

static void long_double_fresnelcos(char **args, npy_intp *dimensions,
                                   npy_intp* steps, void* data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements,
                          Fresnel_Cosine_Taylor_to_Asymptotic_Long_Double);
}

#endif