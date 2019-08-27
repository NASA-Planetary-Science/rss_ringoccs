#ifndef RSS_RINGOCCS_FRESNEL_DIFFRACTION_WRAPPERS_H
#define RSS_RINGOCCS_FRESNEL_DIFFRACTION_WRAPPERS_H

/* cosine and sine are defined here. */
#include <math.h>

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>

/*  Diffraction modeling functions, using Fresnel approximation, found here.  */
#include "__fresnel_diffraction.h"

/*  Various header files required for the C-Python API to work.               */
#include <numpy/ndarraytypes.h>

/*************Square Well Diffraction Using Fresnel Approximation**************/

/******************************************************************************
 *  Function:                                                                 *
 *      complex_float_square_well                                             *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      square well, assuming the Fresnel approximation is valid.             *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) This is a wrapper for Square_Well_Diffraction_Solution_Float,     *
 *          which is defined in __fresnel_diffraction.h. This allows Python   *
 *          to that function, and allows for numpy arrays to be passed in. Ti *
 *          relies on the Numpy UFUNC API, as well as the C-Python API.       *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void complex_float_square_well(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    float *x  =  (float *)args[0];
    float  a  = *(float *)args[1];
    float  b  = *(float *)args[2];
    float  F  = *(float *)args[3];

    /* The output is a pointer to a complex float.                            */
    complex float *out = (complex float *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Float(x[i], a, b, F);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      complex_double_square_well                                            *
 *  Purpose:                                                                  *
 *      Same as complex_float_square_well but for doubles instead of floats.  *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) See complex_float_square_well above for more documentation.       *
 ******************************************************************************/
static void complex_double_square_well(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    double *x  =  (double *)args[0];
    double  a  = *(double *)args[1];
    double  b  = *(double *)args[2];
    double  F  = *(double *)args[3];

    /* The output is a pointer to a complex double.                           */
    complex double *out = (complex double *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      complex_long_double_square_well                                       *
 *  Purpose:                                                                  *
 *      Same as complex_float_square_well but for long doubles.               *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) See complex_float_square_well above for more documentation.       *
 ******************************************************************************/
static void complex_long_double_square_well(char **args, npy_intp *dimensions,
                                            npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    long i;
    long n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    long double *x  =  (long double *)args[0];
    long double  a  = *(long double *)args[1];
    long double  b  = *(long double *)args[2];
    long double  F  = *(long double *)args[3];

    /* The output is a pointer to a complex double.                           */
    complex long double *out = (complex long double *)args[4];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Solution_Long_Double(x[i], a, b, F);
    }
}

/*--------Inverted Square Well Diffraction Using Fresnel Approximation--------*/

static void complex_float_inv_square_well(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];
    float *x  =  (float *)args[0];
    float a   = *(float *)args[1];
    float b   = *(float *)args[2];
    float F   = *(float *)args[3];

    complex float *out = (complex float *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Float(x[i], a, b, F);
    }
}

static void complex_double_inv_square_well(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    double *x  =  (double *)args[0];
    double a   = *(double *)args[1];
    double b   = *(double *)args[2];
    double F   = *(double *)args[3];
    complex double *out = (complex double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

static void complex_long_double_inv_square_well(char **args,
                                                npy_intp *dimensions,
                                                npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    long double *x  =  (long double *)args[0];
    long double a   = *(long double *)args[1];
    long double b   = *(long double *)args[2];
    long double F   = *(long double *)args[3];
    complex long double *out = (complex long double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Inverted_Square_Well_Diffraction_Solution_Double(x[i], a, b, F);
    }
}

/*-------------Phase from Square Well Using Fresnel Approximation-------------*/

static void float_square_well_phase(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    float *x   =  (float *)args[0];
    float a    = *(float *)args[1];
    float b    = *(float *)args[2];
    float F    = *(float *)args[3];
    float *out =  (float *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Float(x[i], a, b, F);
    }
}

static void double_square_well_phase(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    double *x   =  (double *)args[0];
    double a    = *(double *)args[1];
    double b    = *(double *)args[2];
    double F    = *(double *)args[3];
    double *out =  (double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Double(x[i], a, b, F);
    }
}

static void long_double_square_well_phase(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){
    long i;
    long n = dimensions[0];

    long double *x   =  (long double *)args[0];
    long double a    = *(long double *)args[1];
    long double b    = *(long double *)args[2];
    long double F    = *(long double *)args[3];
    long double *out =  (long double *)args[4];

    for (i = 0; i < n; i++) {
        out[i] = Square_Well_Diffraction_Phase_Double(x[i], a, b, F);
    }
}

#endif