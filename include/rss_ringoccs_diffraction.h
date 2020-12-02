/*  Include guard to avoid importing this file twice.                         */
#ifndef _RSS_RINGOCCS_DIFFRACTION_H_
#define _RSS_RINGOCCS_DIFFRACTION_H_

#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  The ringlet and gap modeling functions.                                   */
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Ringlet_Diffraction(double x, double a, double b, double F);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Gap_Diffraction(double x, double a, double b, double F);

/*  Functions for computing the phase of a ringlet.                           */

extern float
rssringoccs_Float_Ringlet_Diffraction_Phase(float x, float a, float b, float F);

extern double
rssringoccs_Double_Ringlet_Diffraction_Phase(double x, double a,
                                             double b, double F);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Square_Wave_Diffraction(double x, double W,
                                            double F, unsigned int N);

extern long double
rssringoccs_LDouble_Ringlet_Diffraction_Phase(long double x, long double a,
                                                 long double b, long double F);

/*  Left and right straightedge modeling tools.                               */
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Right_Straightedge_Diffraction(double x, double edge,
                                                   double F);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Left_Straightedge_Diffraction(double x, double edge,
                                                  double F);

/******************************************************************************
 *--------------------Single Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
rssringoccs_Float_One_Slit_Fraunhofer_Diffraction(float x, float z, float a);

extern double
rssringoccs_Double_One_Slit_Fraunhofer_Diffraction(double x, double z,
                                                   double a);

extern long double
rssringoccs_LDouble_One_Slit_Fraunhofer_Diffraction(long double x,
                                                       long double z,
                                                       long double a);

/******************************************************************************
 *--------------------Double Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
rssringoccs_Float_Two_Slit_Fraunhofer_Diffraction(float x, float z, float a,
                                                  float d, float lambda);

extern double
rssringoccs_Double_Two_Slit_Fraunhofer_Diffraction(double x, double z, double a,
                                                   double d, double lambda);

extern long double
rssringoccs_LDouble_Two_Slit_Fraunhofer_Diffraction(long double x,
                                                       long double z,
                                                       long double a,
                                                       long double d,
                                                       long double lambda);

#endif
