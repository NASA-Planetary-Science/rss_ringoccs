#ifndef RSS_RINGOCCS_PHYSICS_FUNCTIONS_WRAPPERS_H
#define RSS_RINGOCCS_PHYSICS_FUNCTIONS_WRAPPERS_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

/*  Location of the C functions to be wrapped.                                */
#include "__physics_functions.h"

/*  Wrappers for the frequency_to_wavelength function.                        */

static void float_frequency_to_wavelength(char **args, npy_intp *dimensions,
                                          npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Frequency_To_Wavelength_Float(in[i]);
    }
}

static void double_frequency_to_wavelength(char **args, npy_intp *dimensions,
                                           npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Frequency_To_Wavelength_Double(in[i]);
    }
}

static void long_double_frequency_to_wavelength(char **args,
                                                npy_intp *dimensions,
                                                npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Frequency_To_Wavelength_Long_Double(in[i]);
    }
}

/*  Wrappers for the wavelength_to_wavenumber function.                       */

static void float_wavelength_to_wavenumber(char **args, npy_intp *dimensions,
                                           npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Wavelength_To_Wavenumber_Float(in[i]);
    }
}

static void double_wavelength_to_wavenumber(char **args, npy_intp *dimensions,
                                            npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    double *in  = (double *)args[0];
    double *out = (double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Wavelength_To_Wavenumber_Double(in[i]);
    }
}

static void long_double_wavelength_to_wavenumber(char **args,
                                                 npy_intp *dimensions,
                                                 npy_intp *steps, void *data){
    long i;
    long n = dimensions[0];

    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Wavelength_To_Wavenumber_Long_Double(in[i]);
    }
}

/*  Wrappers for the Fresnel scale function.                                  */

static void float_fresnel_scale(char **args, npy_intp *dimensions,
                                npy_intp *steps, void *data){
    npy_intp i, j, k, l, m, dim;

    float *lambda  = (float *)args[0];
    float *d       = (float *)args[1];
    float *phi     = (float *)args[2];
    float *b       = (float *)args[3];
    float *out     = (float *)args[4];

    npy_intp l_step = (steps[0] != 0);
    npy_intp d_step = (steps[1] != 0);
    npy_intp p_step = (steps[2] != 0);
    npy_intp b_step = (steps[3] != 0);

    dim = dimensions[0];

    j = 0;
    k = 0;
    l = 0;
    m = 0;

    for (i = 0; i < dim; i++) {
        out[i] = Fresnel_Scale_Float(lambda[j], d[k], phi[l], b[m]);

        j += l_step;
        k += d_step;
        l += p_step;
        m += b_step;
    }
}

static void double_fresnel_scale(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data){
    npy_intp i, j, k, l, m, dim;

    double *lambda  = (double *)args[0];
    double *d       = (double *)args[1];
    double *phi     = (double *)args[2];
    double *b       = (double *)args[3];
    double *out     = (double *)args[4];

    npy_intp l_step = (steps[0] != 0);
    npy_intp d_step = (steps[1] != 0);
    npy_intp p_step = (steps[2] != 0);
    npy_intp b_step = (steps[3] != 0);

    dim = dimensions[0];


    j = 0;
    k = 0;
    l = 0;
    m = 0;

    for (i = 0; i < dim; i++) {
        out[i] = Fresnel_Scale_Double(lambda[j], d[k], phi[l], b[m]);

        j += l_step;
        k += d_step;
        l += p_step;
        m += b_step;
    }
}

static void long_double_fresnel_scale(char **args, npy_intp *dimensions,
                                      npy_intp *steps, void *data){
    npy_intp i, j, k, l, m, dim;

    long double *lambda  = (long double *)args[0];
    long double *d       = (long double *)args[1];
    long double *phi     = (long double *)args[2];
    long double *b       = (long double *)args[3];
    long double *out     = (long double *)args[4];

    npy_intp l_step = (steps[0] != 0);
    npy_intp d_step = (steps[1] != 0);
    npy_intp p_step = (steps[2] != 0);
    npy_intp b_step = (steps[3] != 0);

    dim = dimensions[0];
    

    j = 0;
    k = 0;
    l = 0;
    m = 0;

    for (i = 0; i < dim; i++) {
        out[i] = Fresnel_Scale_Long_Double(lambda[j], d[k], phi[l], b[m]);

        j += l_step;
        k += d_step;
        l += p_step;
        m += b_step;
    }
}


#endif