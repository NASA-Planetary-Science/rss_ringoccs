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
    long n  = dimensions[0];
    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Frequency_To_Wavelength_Float(in[i]);
    }
}

static void double_frequency_to_wavelength(char **args, npy_intp *dimensions,
                                           npy_intp *steps, void *data){
    long i;
    long n  = dimensions[0];
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
    long n  = dimensions[0];
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
    long n  = dimensions[0];
    float *in  = (float *)args[0];
    float *out = (float *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Wavelength_To_Wavenumber_Float(in[i]);
    }
}

static void double_wavelength_to_wavenumber(char **args, npy_intp *dimensions,
                                            npy_intp *steps, void *data){
    long i;
    long n  = dimensions[0];
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
    long n  = dimensions[0];
    long double *in  = (long double *)args[0];
    long double *out = (long double *)args[1];

    for (i = 0; i < n; i++) {
        out[i] = Wavelength_To_Wavenumber_Long_Double(in[i]);
    }
}

#endif