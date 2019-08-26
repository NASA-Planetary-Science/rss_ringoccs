#ifndef RSS_RINGOCCS_FRESNEL_INTEGRALS_WRAPPER_H
#define RSS_RINGOCCS_FRESNEL_INTEGRALS_WRAPPER_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

#include "__fresnel_integrals.h"

static void double_fresnelsin(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n  = dimensions[0];
    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Fresnel_Sine_Taylor_to_Asymptotic_Func(in[i]);
    }
}

static void double_fresnelcos(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n  = dimensions[0];
    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        /*  The function Fresnel_Cosine_Taylor_to_Asymptotic_Func is defined  *
         *  in _fresnel_cosine.h. Make sure this is in the current directory! */
        out[i] = Fresnel_Cosine_Taylor_to_Asymptotic_Func(in[i]);
    }
}

#endif