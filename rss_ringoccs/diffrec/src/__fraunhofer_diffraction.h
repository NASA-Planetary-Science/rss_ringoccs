#ifndef RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_H
#define RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_H

#include <math.h>
#include "__sinc.h"
#include "__math_constants.h"

/******************************************************************************
 *--------------------Single Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float Single_Slit_Fraunhofer_Diffraction_Float(float x, float z,
                                                      float a);

extern double Single_Slit_Fraunhofer_Diffraction_Double(double x, double z,
                                                        double a);

extern long double Single_Slit_Fraunhofer_Diffraction_Long_Double(
    long double x, long double z, long double a
);

/******************************************************************************
 *--------------------Double Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float Double_Slit_Fraunhofer_Diffraction_Float(float x, float z,
                                                      float a, float d);

extern double Double_Slit_Fraunhofer_Diffraction_Double(double x, double z,
                                                        double a, double d);

extern long double Double_Slit_Fraunhofer_Diffraction_Long_Double(
    long double x, long double z, long double a, long double d
);

#endif