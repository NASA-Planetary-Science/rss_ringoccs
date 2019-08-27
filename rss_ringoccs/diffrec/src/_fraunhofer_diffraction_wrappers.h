#ifndef RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_WRAPPERS_H
#define RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_WRAPPERS_H

#include "__fraunhofer_diffraction.h"

/*  Various header files required for the C-Python API to work.               */
#include <numpy/ndarraytypes.h>

/************Single Slit Diffraction Using Fraunhofer Approximation************/

/******************************************************************************
 *  Function:                                                                 *
 *      float_single_slit_diffraction                                         *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      single slit using the Fraunhofer approximation.                       *
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
 *      1.) This is a wrapper for Single_Slit_Fraunhofer_Diffraction_Float,   *
 *          which is defined in __fraunhofer_diffraction.h. This allows       *
 *          Python to use that function, and allows for numpy arrays to be    *
 *          passed in. Relies on the Numpy UFUNC API and the C-Python API.    *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void float_single_slit_diffraction(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    float *x  =  (float *)args[0];
    float  z  = *(float *)args[1];
    float  a  = *(float *)args[2];

    /* The output is a pointer to a complex float.                            */
    float *out = (float *)args[3];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Single_Slit_Fraunhofer_Diffraction_Float(x[i], z, a);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      double_single_slit_diffraction                                        *
 *  Purpose:                                                                  *
 *      Same as float_single_slit_diffraction for double precision.           *
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
 *      1.) This is a wrapper for Single_Slit_Fraunhofer_Diffraction_Double,  *
 *          which is defined in __fraunhofer_diffraction.h. This allows       *
 *          Python to use that function, and allows for numpy arrays to be    *
 *          passed in. Relies on the Numpy UFUNC API and the C-Python API.    *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void double_single_slit_diffraction(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    double *x  =  (double *)args[0];
    double  z  = *(double *)args[1];
    double  a  = *(double *)args[2];

    /* The output is a pointer to a complex float.                            */
    double *out = (double *)args[3];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Single_Slit_Fraunhofer_Diffraction_Double(x[i], z, a);
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      long_double_single_slit_diffraction                                   *
 *  Purpose:                                                                  *
 *      Same as float_single_slit_diffraction for long double precision.      *
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
 *      1.) This is for Single_Slit_Fraunhofer_Diffraction_Long_Double,       *
 *          which is defined in __fraunhofer_diffraction.h. This allows       *
 *          Python to use that function, and allows for numpy arrays to be    *
 *          passed in. Relies on the Numpy UFUNC API and the C-Python API.    *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void long_double_single_slit_diffraction(char **args,
                                                npy_intp *dimensions,
                                                npy_intp* steps, void* data){

    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    long double *x  =  (long double *)args[0];
    long double  z  = *(long double *)args[1];
    long double  a  = *(long double *)args[2];

    /* The output is a pointer to a complex float.                            */
    long double *out = (long double *)args[3];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Single_Slit_Fraunhofer_Diffraction_Long_Double(x[i], z, a);
    }
}

#endif