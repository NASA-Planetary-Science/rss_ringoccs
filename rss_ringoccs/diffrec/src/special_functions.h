#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>
#include <stdlib.h>

/*  Most functions return floats, not int or int-like data types. To allow    *
 *  most python/numpy data-types functions are created that intake the given  *
 *  data type, convert to double, and then return a double. This preprocessor *
 *  function eliminates excessive and repetitive code. This also ensures the  *
 *  naming convention of functions is consistent throughout the code. For     *
 *  example, consider the Sinc function defined by:                           *
 *                     --                                                     *
 *                    |     sin(x)/x ,      x =/= 0                           *
 *      sinc(x)  =    |                                                       *
 *                    |     1 ,             x = 0                             *
 *                     --                                                     *
 *  This routine is outlined in the sinc file where Sinc_Float, Sinc_Double,  *
 *  and Sinc_Long_Double are defined and use the C Standard Library functions *
 *  sinf, sin, and sinl which compute the sine function at float, double, and *
 *  long double precision accuracy, respectively. If a user at the Python     *
 *  level passes an array of integers it would be convenient to just have a   *
 *  function Sinc_Int that takes in a int, converts to double, and then       *
 *  returns as a double. Similarly if a user passes a char, long, short, etc. *
 *  This preprocessor routine generates all of these functions for a given    *
 *  function FuncName where FuncName_Double has ALREADY been defined in some  *
 *  external file, for example sinc.c.                                        */
#ifdef RSSRINGOCCSNonFloatInputForFloatOutput
#undef RSSRINGOCCSNonFloatInputForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputForFloatOutput(FuncName)\
double FuncName##_Char(char x)          {return FuncName##_Double((double)x);}\
double FuncName##_UChar(unsigned char x){return FuncName##_Double((double)x);}\
double FuncName##_Short(short x)        {return FuncName##_Double((double)x);}\
double FuncName##_UShort(unsigned short x){\
    return FuncName##_Double((double)x);\
}\
double FuncName##_Int(int x)            {return FuncName##_Double((double)x);}\
double FuncName##_UInt(unsigned int x)  {return FuncName##_Double((double)x);}\
double FuncName##_Long(long x)          {return FuncName##_Double((double)x);}\
double FuncName##_ULong(unsigned long x){return FuncName##_Double((double)x);}\
double FuncName##_Long_Long(long long x){return FuncName##_Double((double)x);}\
double FuncName##_ULong_Long(unsigned long long x){\
    return FuncName##_Double((double)x);\
}

/*  This code is for generating the code for the various window functions.    */
#ifdef RSSRINGOCCSNonFloatInputTwoVarForFloatOutput
#undef RSSRINGOCCSNonFloatInputTwoVarForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputTwoVarForFloatOutput(FuncName)\
double FuncName##_Char(char x, double y)                                    \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UChar(unsigned char x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Short(short x, double y)                                  \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UShort(unsigned short x, double y)                        \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Int(int x, double y)                                      \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UInt(unsigned int x, double y)                            \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Long(long x, double y)                                    \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_ULong(unsigned long x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Long_Long(long long x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_ULong_Long(unsigned long long x, double y)                \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}

/*  This code is for generating the code for the straightedge routines.       */
#ifdef RSSRINGOCCSNonFloatInputThreeVarForFloatOutput
#undef RSSRINGOCCSNonFloatInputThreeVarForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputThreeVarForFloatOutput(FuncName, type)      \
type FuncName##_Char(char x, double a, double F)                            \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_UChar(unsigned char x, double a, double F)                  \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_Short(short x, double a, double F)                          \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_UShort(unsigned short x, double a, double F)                \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_Int(int x, double a, double F)                              \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_UInt(unsigned int x, double a, double F)                    \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_Long(long x, double a, double F)                            \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_ULong(unsigned long x, double a, double F)                  \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_Long_Long(long long x, double a, double F)                  \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}                                                                           \
type FuncName##_ULong_Long(unsigned long long x, double a, double F)        \
{                                                                           \
    return FuncName##_Double((double)x, a, F);                              \
}

/*  This code is for generating the code for diffraction modeling functions.  */
#ifdef RSSRINGOCCSNonFloatInputFourVarForFloatOutput
#undef RSSRINGOCCSNonFloatInputFourVarForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputFourVarForFloatOutput(FuncName, type)\
type FuncName##_Char(char x, double a, double b, double F)                  \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_UChar(unsigned char x, double a, double b, double F)        \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_Short(short x, double a, double b, double F)                \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_UShort(unsigned short x, double a, double b, double F)      \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_Int(int x, double a, double b, double F)                    \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_UInt(unsigned int x, double a, double b, double F)          \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_Long(long x, double a, double b, double F)                  \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_ULong(unsigned long x, double a, double b, double F)        \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type FuncName##_Long_Long(long long x, double a, double b, double F)        \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}                                                                           \
type                                                                        \
FuncName##_ULong_Long(unsigned long long x, double a, double b, double F)   \
{                                                                           \
    return FuncName##_Double((double)x, a, b, F);                           \
}

/*  These lines are repeated over and over in this header file to define the  *
 *  various math functions. Use this preprocessor function to save lines of   *
 *  code, stay consistent with naming conventions, and add clarity.           */
#ifdef RSSRINGOCCSGenerateExternFunctions
#undef RSSRINGOCCSGenerateExternFunctions
#endif

/*  Macros to avoid repetition and ensure consistency in naming convention.   */
#ifdef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSTwoVarWindowFuncExtern
#endif

#ifdef RSSRINGOCCSDiffractionModelingExtern
#undef RSSRINGOCCSDiffractionModelingExtern
#endif

#ifdef RSSRINGOCCSStraightedgeModelingExtern
#undef RSSRINGOCCSStraightedgeModelingExtern
#endif

#define RSSRINGOCCSGenerateExternFunctions(FuncName)                        \
extern float        FuncName##_Float(float x);                              \
extern double       FuncName##_Double(double x);                            \
extern long double  FuncName##_Long_Double(long double x);                  \
extern double       FuncName##_Char(char x);                                \
extern double       FuncName##_UChar(unsigned char x);                      \
extern double       FuncName##_Short(short x);                              \
extern double       FuncName##_UShort(unsigned short x);                    \
extern double       FuncName##_Int(int x);                                  \
extern double       FuncName##_UInt(unsigned int x);                        \
extern double       FuncName##_Long(long x);                                \
extern double       FuncName##_ULong(unsigned long x);                      \
extern double       FuncName##_Long_Long(long long x);                      \
extern double       FuncName##_ULong_Long(unsigned long long x);

#define RSSRINGOCCSTwoVarWindowFuncExtern(FuncName)                         \
extern float        FuncName##_Float(float x, float W);                     \
extern double       FuncName##_Double(double x, double W);                  \
extern long double  FuncName##_Long_Double(long double x, long double W);   \
extern double       FuncName##_Char(char x, double W);                      \
extern double       FuncName##_UChar(unsigned char x, double W);            \
extern double       FuncName##_Short(short x, double W);                    \
extern double       FuncName##_UShort(unsigned short x, double W);          \
extern double       FuncName##_Int(int x, double W);                        \
extern double       FuncName##_UInt(unsigned int x, double W);              \
extern double       FuncName##_Long(long x, double W);                      \
extern double       FuncName##_ULong(unsigned long x, double W);            \
extern double       FuncName##_Long_Long(long long x, double W);            \
extern double       FuncName##_ULong_Long(unsigned long long x, double W);

#define RSSRINGOCCSThreeVarWindowFuncExtern(FuncName)                       \
extern float  FuncName##_Float(float x, float W, float alpha);              \
extern double FuncName##_Double(double x, double W, double alpha);          \
extern long double                                                          \
FuncName##_Long_Double(long double x, long double W, long double alpha);    \
extern double FuncName##_Char(char x, double W, double alpha);              \
extern double FuncName##_UChar(unsigned char x, double W, double alpha);    \
extern double FuncName##_Short(short x, double W, double alpha);            \
extern double FuncName##_UShort(unsigned short x, double W, double alpha);  \
extern double FuncName##_Int(int x, double W, double alpha);                \
extern double FuncName##_UInt(unsigned int x, double W, double alpha);      \
extern double FuncName##_Long(long x, double W, double alpha);              \
extern double FuncName##_ULong(unsigned long x, double W, double alpha);    \
extern double FuncName##_Long_Long(long long x, double W, double alpha);    \
extern double                                                               \
FuncName##_ULong_Long(unsigned long long x, double W, double alpha);

#define RSSRINGOCCSDiffractionModelingExtern(FuncName, type)                \
type float                                                                  \
FuncName##_Float(float x, float a, float b, float F);                       \
type double                                                                 \
FuncName##_Double(double x, double a, double b, double F);                  \
type long double                                                            \
FuncName##_Long_Double(long double x, long double a,                        \
                       long double b, long double F);                       \
type double                                                                 \
FuncName##_Char(char x, double a, double b, double F);                      \
type double                                                                 \
FuncName##_UChar(unsigned char x, double a, double b, double F);            \
type double                                                                 \
FuncName##_Short(short x, double a, double b, double F);                    \
type double                                                                 \
FuncName##_UShort(unsigned short x, double a, double b, double F);          \
type double                                                                 \
FuncName##_Int(int x, double a, double b, double F);                        \
type double                                                                 \
FuncName##_UInt(unsigned int x, double a, double b, double F);              \
type double                                                                 \
FuncName##_Long(long x, double a, double b, double F);                      \
type double                                                                 \
FuncName##_ULong(unsigned long x, double a, double b, double F);            \
type double                                                                 \
FuncName##_Long_Long(long long x, double a, double b, double F);            \
type double                                                                 \
FuncName##_ULong_Long(unsigned long long x, double a, double b, double F);

#define RSSRINGOCCSStraightedgeModelingExtern(FuncName)                     \
extern complex float    FuncName##_Float(float x, float edge, float F);     \
extern complex double   FuncName##_Double(double x, double edge, double F); \
extern complex long double                                                  \
FuncName##_Long_Double(long double x, long double edge, long double F);     \
extern complex double   FuncName##_Char(char x, double edge, double F);     \
extern complex double                                                       \
FuncName##_UChar(unsigned char x, double edge, double F);                   \
extern complex double   FuncName##_Short(short x, double edge, double F);   \
extern complex double                                                       \
FuncName##_UShort(unsigned short x, double edge, double F);                 \
extern complex double   FuncName##_Int(int x, double edge, double F);       \
extern complex double                                                       \
FuncName##_UInt(unsigned int x, double edge, double F);                     \
extern complex double   FuncName##_Long(long x, double edge, double F);     \
extern complex double                                                       \
FuncName##_ULong(unsigned long x, double edge, double F);                   \
extern complex double                                                       \
FuncName##_Long_Long(long long x, double edge, double F);                   \
extern complex double                                                       \
FuncName##_ULong_Long(unsigned long long x, double edge, double F);

/*  Generate extern function names for all of the math functions.             */
RSSRINGOCCSGenerateExternFunctions(BesselJ0);
RSSRINGOCCSGenerateExternFunctions(BesselI0);
RSSRINGOCCSGenerateExternFunctions(LambertW);
RSSRINGOCCSGenerateExternFunctions(Sinc);
RSSRINGOCCSGenerateExternFunctions(Wavelength_To_Wavenumber);
RSSRINGOCCSGenerateExternFunctions(Frequency_To_Wavelength);
RSSRINGOCCSGenerateExternFunctions(Resolution_Inverse);
RSSRINGOCCSGenerateExternFunctions(Fresnel_Cos);
RSSRINGOCCSGenerateExternFunctions(Fresnel_Sin);

/******************************************************************************
 *------------------------Complex Fresnel Integral----------------------------*
 ******************************************************************************/
extern complex double Fresnel_Taylor_to_Asymptotic_Double(double x);
extern complex double Fresnel_While_to_Asymptotic_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Three_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Four_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Six_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Eight_Func(double x);

/*  Kaiser-Bessel function with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_0);
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_5);
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_3_5);

/*  Modified Kaiser-Bessel with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_0);
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_5);
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_3_5);

/*  Rectangular and Squared Cosine window functions.                          */
RSSRINGOCCSTwoVarWindowFuncExtern(Rect_Window);
RSSRINGOCCSTwoVarWindowFuncExtern(Coss_Window);

/*  Modified and unmodified Kaiser-Bessel function with arbitrary alpha.      */
RSSRINGOCCSThreeVarWindowFuncExtern(Kaiser_Bessel_Al);
RSSRINGOCCSThreeVarWindowFuncExtern(Modified_Kaiser_Bessel_Al);

/*  The ringlet and gap modeling functions.                                   */
RSSRINGOCCSDiffractionModelingExtern(Ringlet_Diffraction, extern complex);
RSSRINGOCCSDiffractionModelingExtern(Gap_Diffraction, extern complex);
RSSRINGOCCSDiffractionModelingExtern(Ringlet_Diffraction_Phase, extern);

/*  Left and right straightedge modeling tools.                               */
RSSRINGOCCSStraightedgeModelingExtern(Right_Straightedge_Diffraction);
RSSRINGOCCSStraightedgeModelingExtern(Left_Straightedge_Diffraction);

#undef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSThreeVarWindowFuncExtern
#undef RSSRINGOCCSDiffractionModelingExtern
#undef RSSRINGOCCSStraightedgeModelingExtern

extern void Legendre_Polynomials(double *legendre_p, double x, int order);

extern void
Alt_Legendre_Polynomials(double *poly, double *legendre_p, int order);

extern void
Fresnel_Kernel_Coefficients(double *fresnel_ker_coeffs, double *legendre_p,
                            double *alt_legendre_p, double Legendre_Coeff,
                            int order);

extern float        Max_Float(float *arr, long n_elements);
extern double       Max_Double(double *arr, long n_elements);
extern long         double Max_Long_Double(long double *arr, long n_elements);
extern short        Max_Short(short *arr, long n_elements);
extern int          Max_Int(int *arr, long n_elements);
extern long         Max_Long(long *arr, long n_elements);
extern long long    Max_Long_Long(long long *arr, long n_elements);

extern float        Min_Float(float *arr, long n_elements);
extern double       Min_Double(double *arr, long n_elements);
extern long         double Min_Long_Double(long double *arr, long n_elements);
extern short        Min_Short(short *arr, long n_elements);
extern int          Min_Int(int *arr, long n_elements);
extern long         Min_Long(long *arr, long n_elements);
extern long long    Min_Long_Long(long long *arr, long n_elements);

extern float        Normeq_Float(float *w_func, long n_elements);
extern double       Normeq_Double(double *w_func, long n_elements);
extern long double  Normeq_Long_Double(long double *w_func, long n_elements);
extern double       Normeq_Short(short *w_func, long n_elements);
extern double       Normeq_Int(int *w_func, long n_elements);
extern double       Normeq_Long(long *w_func, long n_elements);
extern double       Normeq_Long_Long(long long *w_func, long n_elements);

/*  Window Normalization Functions                                            */
extern float
Window_Normalization_Float(float *ker, long dim, float dx, float f_scale);

extern double
Window_Normalization_Double(double *ker, long dim, double dx, double f_scale);

extern long double
Window_Normalization_Long_Double(long double *ker, long dim,
                                 long double dx, long double f_scale);

extern double
Window_Normalization_Short(short *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Int(int *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long(long *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long_Long(long long *ker, long dim,
                               double dx, double f_scale);

extern float
Window_Normalization_Complex_Float(complex float *ker, long dim,
                                   float dx, float f_scale);

extern double
Window_Normalization_Complex_Double(complex double *ker, long dim,
                                    double dx, double f_scale);

extern long double
Window_Normalization_Complex_Long_Double(complex long double *ker, long dim,
                                         long double dx, long double f_scale);

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

/*------------Square Wave Diffraction Using Fresnel Approximation-------------*/
extern complex double
Square_Wave_Diffraction_Double(double x, double W, double F, long N);

#endif
