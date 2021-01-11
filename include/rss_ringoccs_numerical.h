/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

/*  Include guard for this file to prevent including this twice.              */
#ifndef _RSS_RINGOCCS_NUMERICAL_H_
#define _RSS_RINGOCCS_NUMERICAL_H_

#include <rss_ringoccs/include/rss_ringoccs_complex.h>

extern float
rssringoccs_Newton_Raphson_Float(float x, float (*f)(float),
                                 float (*f_prime)(float),
                                 unsigned int max_iters);

extern double
rssringoccs_Newton_Raphson_Double(double x, double (*f)(double),
                                  double (*f_prime)(double),
                                  unsigned int max_iters);

extern long double
rssringoccs_Newton_Raphson_LDouble(long double x,
                                   long double (*f)(long double),
                                   long double (*f_prime)(long double),
                                   unsigned int max_iters);

extern rssringoccs_ComplexDouble
rssringoccs_Newton_Raphson_Complex(
    rssringoccs_ComplexDouble z,
    rssringoccs_ComplexDouble (*f)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_prime)(rssringoccs_ComplexDouble),
    unsigned int max_iters
);

extern float
rssringoccs_Halleys_Method_Float(float x, float (*f)(float),
                                 float (*f_prime)(float),
                                 float (*f_2prime)(float),
                                 unsigned int max_iters);

extern double
rssringoccs_Halleys_Method_Double(double x, double (*f)(double),
                                  double (*f_prime)(double),
                                  double (*f_2prime)(double),
                                  unsigned int max_iters);

extern long double
rssringoccs_Halleys_Method_LDouble(long double x,
                                   long double (*f)(long double),
                                   long double (*f_prime)(long double),
                                   long double (*f_2prime)(long double),
                                   unsigned int max_iters);

extern rssringoccs_ComplexDouble
rssringoccs_Halleys_Method_Complex(
    rssringoccs_ComplexDouble z,
    rssringoccs_ComplexDouble (*f)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_prime)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_2prime)(rssringoccs_ComplexDouble),
    unsigned int max_iters
);

extern float
rssringoccs_Float_Five_Point_Derivative(float (*f)(float), float x, float h);

extern double
rssringoccs_Double_Five_Point_Derivative(double (*f)(double),
                                         double x, double h);

extern long double
rssringoccs_LDouble_Five_Point_Derivative(long double (*f)(long double),
                                             long double x, long double h);


extern rssringoccs_ComplexDouble
rssringoccs_Newton_Raphson_CDouble_Poly_Real(
    rssringoccs_ComplexDouble z, double *coeffs, unsigned int degree,
    unsigned int max_iters
);

#endif
/*  End of include guard: #ifndef _RSS_RINGOCCS_NUMERICAL_H_                  */
