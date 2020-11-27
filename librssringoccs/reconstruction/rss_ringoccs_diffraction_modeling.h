/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_DIFFRACTION_MODELING_H
#define RSS_RINGOCCS_DIFFRACTION_MODELING_H

#include <complex.h>

#ifdef RSSRINGOCCSDiffractionModelingExtern
#undef RSSRINGOCCSDiffractionModelingExtern
#endif

#ifdef RSSRINGOCCSStraightedgeModelingExtern
#undef RSSRINGOCCSStraightedgeModelingExtern
#endif

#define RSSRINGOCCSDiffractionModelingExtern(FuncName, type)                   \
type float                                                                     \
FuncName##_Float(float x, float a, float b, float F);                          \
type double                                                                    \
FuncName##_Double(double x, double a, double b, double F);                     \
type long double                                                               \
FuncName##_Long_Double(long double x, long double a,                           \
                       long double b, long double F);

#define RSSRINGOCCSStraightedgeModelingExtern(FuncName)                        \
extern complex float    FuncName##_Float(float x, float edge, float F);        \
extern complex double   FuncName##_Double(double x, double edge, double F);    \
extern complex long double                                                     \
FuncName##_Long_Double(long double x, long double edge, long double F);

/*  The ringlet and gap modeling functions.                                   */
RSSRINGOCCSDiffractionModelingExtern(Ringlet_Diffraction, extern complex);
RSSRINGOCCSDiffractionModelingExtern(Gap_Diffraction, extern complex);
RSSRINGOCCSDiffractionModelingExtern(Ringlet_Diffraction_Phase, extern);

/*  Left and right straightedge modeling tools.                               */
RSSRINGOCCSStraightedgeModelingExtern(Right_Straightedge_Diffraction);
RSSRINGOCCSStraightedgeModelingExtern(Left_Straightedge_Diffraction);

#undef RSSRINGOCCSDiffractionModelingExtern
#undef RSSRINGOCCSStraightedgeModelingExtern

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

#endif