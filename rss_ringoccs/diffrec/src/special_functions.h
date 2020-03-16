#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>
#include <stdlib.h>

/*  compute_norm_eq, max and min found here. math.h included here as well.    */
#include "__fresnel_diffraction.h"
#include "__math_functions.h"
#include "__get_array.h"
#include "__window_functions.h"
#include "__diffraction_functions.h"

/*  Make sure these names are available.                                      */
#ifdef __get_one_real_from_one_real
#undef __get_one_real_from_one_real
#endif

#ifdef __get_one_real_from_two_real
#undef __get_one_real_from_two_real
#endif

#ifdef __get_one_real_from_three_real
#undef __get_one_real_from_three_real
#endif

#ifdef __get_complex_from_four_real
#undef __get_complex_from_four_real
#endif

#ifdef __get_complex_from_three_real
#undef __get_complex_from_three_real
#endif

/*  To avoid repeating the same code over and over again, define these macros *
 *  to be used for looping over functions.                                    */

#define __get_one_real_from_one_real(x, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i]);\
    }\
})

#define __get_one_real_from_two_real(x1, x2, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x1[i], x2);\
    }\
})

#define __get_one_real_from_three_real(x1, x2, x3, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x1[i], x2, x3);\
    }\
})

#define __get_complex_from_three_real(x, a, F, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i], a, F);\
    }\
})

#define __get_complex_from_four_real(x, a, b, F, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i], a, b, F);\
    }\
})

/******************************************************************************
 *--------------------Single Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
Single_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a);

extern double
Single_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a);

extern long double
Single_Slit_Fraunhofer_Diffraction_Long_Double(long double x, long double z,
                                               long double a);

/******************************************************************************
 *--------------------Double Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
Double_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a,
                                         float d, float lambda);

extern double
Double_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a,
                                          double d, double lambda);

extern long double
Double_Slit_Fraunhofer_Diffraction_Long_Double(long double x, long double z,
                                               long double a, long double d,
                                               long double lambda);

/******************************************************************************
 *------------------------------Fresnel Scale---------------------------------*
 ******************************************************************************/

extern float
Fresnel_Scale_Float(float lambda, float d, float phi, float b);

extern double
Fresnel_Scale_Double(double lambda, double d, double phi, double b);

extern long double
Fresnel_Scale_Long_Double(long double lambda, long double d,
                          long double phi, long double b);

/******************************************************************************
 *--------------------------Frequency to Wavelength---------------------------*
 ******************************************************************************/

extern float        Frequency_To_Wavelength_Float(float frequency);
extern double       Frequency_To_Wavelength_Double(double frequency);
extern long double  Frequency_To_Wavelength_Long_Double(long double frequency);

/******************************************************************************
 *-------------------------Wavelength to Wavenumber---------------------------*
 ******************************************************************************/

extern float        Wavelength_To_Wavenumber_Float(float frequency);
extern double       Wavelength_To_Wavenumber_Double(double frequency);
extern long double  Wavelength_To_Wavenumber_Long_Double(long double frequency);

/******************************************************************************
 *-----------------------------------Where------------------------------------*
 ******************************************************************************/

extern long **
Where_Greater_Char(char *data, long dim, double threshold);

extern long **
Where_Greater_UChar(unsigned char *data, long dim, double threshold);

extern long **
Where_Greater_Short(short *data, long dim, double threshold);

extern long **
Where_Greater_UShort(unsigned short *data, long dim, double threshold);

extern long **
Where_Greater_Int(int *data, long dim, double threshold);

extern long **
Where_Greater_UInt(unsigned int *data, long dim, double threshold);

extern long **
Where_Greater_Long(long *data, long dim, double threshold);

extern long **
Where_Greater_ULong(unsigned long *data, long dim, double threshold);

extern long **
Where_Greater_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Greater_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **
Where_Greater_Float(float *data, long dim, double threshold);

extern long **
Where_Greater_Double(double *data, long dim, double threshold);

extern long **
Where_Greater_Long_Double(long double *data, long dim, long double threshold);

extern long **
Where_Lesser_Char(char *data, long dim, double threshold);

extern long **
Where_Lesser_UChar(unsigned char *data, long dim, double threshold);

extern long **
Where_Lesser_Short(short *data, long dim, double threshold);

extern long **
Where_Lesser_UShort(unsigned short *data, long dim, double threshold);

extern long **
Where_Lesser_Int(int *data, long dim, double threshold);

extern long **
Where_Lesser_UInt(unsigned int *data, long dim, double threshold);

extern long **
Where_Lesser_Long(long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong(unsigned long *data, long dim, double threshold);

extern long **
Where_Lesser_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **
Where_Lesser_Float(float *data, long dim, double threshold);

extern long **
Where_Lesser_Double(double *data, long dim, double threshold);

extern long **
Where_Lesser_Long_Double(long double *data, long dim, long double threshold);

#endif