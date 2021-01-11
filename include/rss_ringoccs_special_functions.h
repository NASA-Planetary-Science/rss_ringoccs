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

/*  typedef for the window function pointers. Window functions take in two    *
 *  doubles (the x-value and the window width), and return a double.          *
 *  The function pointer works as follows:                                    *
 *      return_type  (*type_name)(type_var1, type_var2, ...)                  *
 *  So, let's typedef this for the window function.                           */
typedef double (*rssringoccs_window_func)(double, double);

/*  As a side comment, the FresT function pointer takes a different number of *
 *  variables depending on which method of diffraction correction is being    *
 *  performed, so we can't just typedef it here. We'll need to declare it     *
 *  individually for each diffraction correction method instead.              */

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
extern float rssringoccs_Float_##FuncName(float x);                            \
extern double rssringoccs_Double_##FuncName(double x);                         \
extern long double  rssringoccs_LDouble_##FuncName(long double x);             \
extern rssringoccs_ComplexDouble                                               \
rssringoccs_CDouble_##FuncName(rssringoccs_ComplexDouble x);

#define RSSRINGOCCSTwoVarWindowFuncExtern(FuncName)                            \
extern float rssringoccs_Float_##FuncName(float x, float W);                   \
extern double rssringoccs_Double_##FuncName(double x, double W);               \
extern long double                                                             \
rssringoccs_LDouble_##FuncName(long double x, long double W);

#define RSSRINGOCCSThreeVarWindowFuncExtern(FuncName)                          \
extern float  rssringoccs_Float_##FuncName(float x, float W, float alpha);     \
extern double rssringoccs_Double_##FuncName(double x, double W, double alpha); \
extern long double                                                             \
rssringoccs_LDouble_##FuncName(long double x, long double W,                \
                                  long double alpha);

/*  Generate extern function names for all of the math functions.             */
RSSRINGOCCSGenerateExternFunctions(Bessel_J0)
RSSRINGOCCSGenerateExternFunctions(Bessel_I0)
RSSRINGOCCSGenerateExternFunctions(LambertW)
RSSRINGOCCSGenerateExternFunctions(Wavelength_To_Wavenumber)
RSSRINGOCCSGenerateExternFunctions(Frequency_To_Wavelength)
RSSRINGOCCSGenerateExternFunctions(Fresnel_Cos)
RSSRINGOCCSGenerateExternFunctions(Fresnel_Sin)

extern float rssringoccs_Float_Resolution_Inverse(float x);
extern double rssringoccs_Double_Resolution_Inverse(double x);
extern long double rssringoccs_LDouble_Resolution_Inverse(long double x);

extern rssringoccs_ComplexDouble rssringoccs_Complex_Fresnel_Integral(double x);

/*  Kaiser-Bessel function with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_0)
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_5)
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_3_5)

/*  Modified Kaiser-Bessel with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_0)
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_5)
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_3_5)

/*  Rectangular and Squared Cosine window functions.                          */
RSSRINGOCCSTwoVarWindowFuncExtern(Rect_Window)
RSSRINGOCCSTwoVarWindowFuncExtern(Coss_Window)

/*  Modified and unmodified Kaiser-Bessel function with arbitrary alpha.      */
RSSRINGOCCSThreeVarWindowFuncExtern(Kaiser_Bessel)
RSSRINGOCCSThreeVarWindowFuncExtern(Modified_Kaiser_Bessel)

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
Fresnel_Scale_LDouble(long double lambda, long double d,
                         long double phi, long double b);

extern float
rssringoccs_Max_Float(float *arr, long n_elements);

extern double
rssringoccs_Max_Double(double *arr, long n_elements);

extern long double
rssringoccs_Max_LDouble(long double *arr, long n_elements);

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
rssringoccs_Min_LDouble(long double *arr, long n_elements);

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
rssringoccs_Normeq_LDouble(long double *w_func, long n_elements);

extern double
rssringoccs_Normeq_Short(short *w_func, long n_elements);

extern double
rssringoccs_Normeq_Int(int *w_func, long n_elements);

extern double
rssringoccs_Normeq_Long(long *w_func, long n_elements);


/*  Window Normalization Functions                                            */
extern float
rssringoccs_Float_Window_Normalization(float *ker, long dim,
                                       float dx, float f_scale);

extern double
rssringoccs_Double_Window_Normalization(double *ker, long dim,
                                        double dx, double f_scale);

extern long double
rssringoccs_LDouble_Window_Normalization(long double *ker,
                                            long dim, long double dx,
                                            long double f_scale);

extern double
rssringoccs_Complex_Window_Normalization(rssringoccs_ComplexDouble *ker,
                                         long dim, double dx, double f_scale);

/******************************************************************************
 *-----------------------------------Where------------------------------------*
 ******************************************************************************/

#define WhereExternFuncs(FuncName)                                             \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Char(char *data, unsigned long dim, double threshold);   \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_UChar(unsigned char *data, unsigned long dim,            \
                                     double threshold);                        \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Short(short *data, unsigned long dim, double threshold); \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_UShort(unsigned short *data, unsigned long dim,          \
                                      double threshold);                       \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Int(int *data, unsigned long dim, double threshold);     \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_UInt(unsigned int *data, unsigned long dim,              \
                                    double threshold);                         \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Long(long *data, unsigned long dim, double threshold);   \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_ULong(unsigned long *data, unsigned long dim,            \
                                     double threshold);                        \
                                                                               \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Float(float *data, unsigned long dim, float threshold);  \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_Double(double *data, unsigned long dim,                  \
                                      double threshold);                       \
                                                                               \
extern unsigned long **                                                                 \
rssringoccs_Where_##FuncName##_LDouble(long double *data, unsigned long dim,         \
                                          long double threshold);

WhereExternFuncs(Lesser)
WhereExternFuncs(Greater)

extern unsigned long **
rssringoccs_Where_LesserGreater_Char(char *data, unsigned long dim,
                                     double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_UChar(unsigned char *data, unsigned long dim,
                                      double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_Short(short *data, unsigned long dim,
                                      double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_UShort(unsigned short *data, unsigned long dim,
                                       double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_Int(int *data, unsigned long dim,
                                    double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_UInt(unsigned int *data, unsigned long dim,
                                     double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_Long(long *data, unsigned long dim,
                                     double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_ULong(unsigned long *data, unsigned long dim,
                                      double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_Float(float *data, unsigned long dim,
                                      float lower, float upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_Double(double *data, unsigned long dim,
                                       double lower, double upper);

extern unsigned long **
rssringoccs_Where_LesserGreater_LDouble(long double *data, unsigned long dim,
                                           long double lower,
                                           long double upper);

#endif
