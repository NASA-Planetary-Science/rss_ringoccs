#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Boolean data type defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

/*  Prototypes for these functions declared here.                             */
#include "rss_ringoccs_special_functions.h"

#define RectNormEQ 1.0
#define CossNormEQ 1.5
#define KB20NormEQ 1.49634231
#define KB25NormEQ 1.65191895
#define KB35NormEQ 1.92844639
#define KBMD20NormEQ 1.52048382
#define KBMD25NormEQ 1.65994438
#define KBMD35NormEQ 1.52048382

/*  These lines are repeated over and over in this header file to define the  *
 *  various math functions. Use this preprocessor function to save lines of   *
 *  code, stay consistent with naming conventions, and add clarity.           */
#ifdef RSSRINGOCCSGenerateExternFunctions
#undef RSSRINGOCCSGenerateExternFunctions
#endif

#ifdef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSTwoVarWindowFuncExtern
#endif

#ifdef RSSRINGOCCSThreeVarWindowFuncExtern
#undef RSSRINGOCCSThreeVarWindowFuncExtern
#endif

#define RSSRINGOCCSGenerateExternFunctions(FuncName)                           \
extern float        FuncName##_Float(float x);                                 \
extern double       FuncName##_Double(double x);                               \
extern long double  FuncName##_LongDouble(long double x);

#define RSSRINGOCCSTwoVarWindowFuncExtern(FuncName)                            \
extern float        FuncName##_Float(float x, float W);                        \
extern double       FuncName##_Double(double x, double W);                     \
extern long double  FuncName##_LongDouble(long double x, long double W);

#define RSSRINGOCCSThreeVarWindowFuncExtern(FuncName)                          \
extern float  FuncName##_Float(float x, float W, float alpha);                  \
extern double FuncName##_Double(double x, double W, double alpha);              \
extern long double                                                             \
FuncName##_LongDouble(long double x, long double W, long double alpha);

/*  Generate extern function names for all of the math functions.             */
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Bessel_J0)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Bessel_I0)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_LambertW)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Sinc)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Wavelength_To_Wavenumber)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Frequency_To_Wavelength)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Resolution_Inverse)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Fresnel_Cos)
RSSRINGOCCSGenerateExternFunctions(rssringoccs_Fresnel_Sin)

extern rssringoccs_ComplexDouble
rssringoccs_Bessel_I0_ComplexDouble(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble rssringoccs_Fresnel_Integral_Double(double x);

/*  Kaiser-Bessel function with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Kaiser_Bessel_2_0)
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Kaiser_Bessel_2_5)
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Kaiser_Bessel_3_5)

/*  Modified Kaiser-Bessel with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Modified_Kaiser_Bessel_2_0)
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Modified_Kaiser_Bessel_2_5)
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Modified_Kaiser_Bessel_3_5)

/*  Rectangular and Squared Cosine window functions.                          */
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Rect_Window)
RSSRINGOCCSTwoVarWindowFuncExtern(rssringoccs_Coss_Window)

/*  Modified and unmodified Kaiser-Bessel function with arbitrary alpha.      */
RSSRINGOCCSThreeVarWindowFuncExtern(rssringoccs_Kaiser_Bessel)
RSSRINGOCCSThreeVarWindowFuncExtern(rssringoccs_Modified_Kaiser_Bessel)

#undef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSThreeVarWindowFuncExtern
#undef RSSRINGOCCSGenerateExternFunctions

extern void
rssringoccs_Legendre_Polynomials(double *legendre_p, double x, int order);

extern void
rssringoccs_Alt_Legendre_Polynomials(double *poly,
                                     double *legendre_p, int order);

extern void
rssringoccs_Fresnel_Kernel_Coefficients(double *fresnel_ker_coeffs,
                                        double *legendre_p,
                                        double *alt_legendre_p,
                                        double Legendre_Coeff, int order);

/******************************************************************************
 *------------------------------Fresnel Scale---------------------------------*
 ******************************************************************************/

extern float
Fresnel_Scale_Float(float lambda, float d, float phi, float b);

extern double
Fresnel_Scale_Double(double lambda, double d, double phi, double b);

extern long double
Fresnel_Scale_LongDouble(long double lambda, long double d,
                         long double phi, long double b);

extern float
rssringoccs_Max_Float(float *arr, long n_elements);

extern double
rssringoccs_Max_Double(double *arr, long n_elements);

extern long double
rssringoccs_Max_LongDouble(long double *arr, long n_elements);

extern char
rssringoccs_Max_Char(char *arr, long n_elements);

extern unsigned char
rssringoccs_Max_UChar(unsigned char *arr, long n_elements);

extern short
rssringoccs_Max_Short(short *arr, long n_elements);

extern unsigned short
rssringoccs_Max_UShort(unsigned short *arr, long n_elements);

extern int
rssringoccs_Max_Int(int *arr, long n_elements);

extern unsigned int
rssringoccs_Max_UInt(unsigned int *arr, long n_elements);

extern long
rssringoccs_Max_Long(long *arr, long n_elements);

extern unsigned long
rssringoccs_Max_ULong(unsigned long *arr, long n_elements);


extern float
rssringoccs_Min_Float(float *arr, long n_elements);

extern double
rssringoccs_Min_Double(double *arr, long n_elements);

extern long double
rssringoccs_Min_LongDouble(long double *arr, long n_elements);

extern char
rssringoccs_Min_Char(char *arr, long n_elements);

extern unsigned char
rssringoccs_Min_UChar(unsigned char *arr, long n_elements);

extern short
rssringoccs_Min_Short(short *arr, long n_elements);

extern unsigned short
rssringoccs_Min_UShort(unsigned short *arr, long n_elements);

extern int
rssringoccs_Min_Int(int *arr, long n_elements);

extern unsigned int
rssringoccs_Min_UInt(unsigned int *arr, long n_elements);

extern long
rssringoccs_Min_Long(long *arr, long n_elements);

extern unsigned long
rssringoccs_Min_ULong(unsigned long *arr, long n_elements);


extern float
rssringoccs_Normeq_Float(float *w_func, long n_elements);

extern double
rssringoccs_Normeq_Double(double *w_func, long n_elements);

extern long double
rssringoccs_Normeq_LongDouble(long double *w_func, long n_elements);

extern double
rssringoccs_Normeq_Short(short *w_func, long n_elements);

extern double
rssringoccs_Normeq_Int(int *w_func, long n_elements);

extern double
rssringoccs_Normeq_Long(long *w_func, long n_elements);


/*  Window Normalization Functions                                            */
extern float
rssringoccs_Window_Normalization_Float(float *ker, long dim,
                                       float dx, float f_scale);

extern double
rssringoccs_Window_Normalization_Double(double *ker, long dim,
                                        double dx, double f_scale);

extern long double
rssringoccs_Window_Normalization_LongDouble(long double *ker,
                                            long dim, long double dx,
                                            long double f_scale);

extern double
rssringoccs_Window_Normalization_ComplexDouble(rssringoccs_ComplexDouble *ker,
                                               long dim, double dx,
                                               double f_scale);

/******************************************************************************
 *-----------------------------------Where------------------------------------*
 ******************************************************************************/

#define WhereExternFuncs(FuncName)                                             \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Char(char *data, long dim, double threshold);   \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_UChar(unsigned char *data, long dim,            \
                                     double threshold);                        \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Short(short *data, long dim, double threshold); \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_UShort(unsigned short *data, long dim,          \
                                      double threshold);                       \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Int(int *data, long dim, double threshold);     \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_UInt(unsigned int *data, long dim,              \
                                    double threshold);                         \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Long(long *data, long dim, double threshold);   \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_ULong(unsigned long *data, long dim,            \
                                     double threshold);                        \
                                                                               \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Float(float *data, long dim, float threshold);  \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_Double(double *data, long dim,                  \
                                      double threshold);                       \
                                                                               \
extern long **                                                                 \
rssringoccs_Where_##FuncName##_LongDouble(long double *data, long dim,         \
                                          long double threshold);

WhereExternFuncs(Lesser)
WhereExternFuncs(Greater)

extern long **
rssringoccs_Where_LesserGreater_Char(char *data, long dim,
                                     double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_UChar(unsigned char *data, long dim,
                                      double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_Short(short *data, long dim,
                                      double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_UShort(unsigned short *data, long dim,
                                       double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_Int(int *data, long dim,
                                    double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_UInt(unsigned int *data, long dim,
                                     double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_Long(long *data, long dim,
                                     double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_ULong(unsigned long *data, long dim,
                                      double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_Float(float *data, long dim,
                                      float lower, float upper);

extern long **
rssringoccs_Where_LesserGreater_Double(double *data, long dim,
                                       double lower, double upper);

extern long **
rssringoccs_Where_LesserGreater_LongDouble(long double *data, long dim,
                                           long double lower,
                                           long double upper);

#endif
