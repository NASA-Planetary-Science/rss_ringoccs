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

#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Include guard for this file to prevent including this twice.              */
#ifndef _RSS_RINGOCCS_NUMERICAL_H_
#define _RSS_RINGOCCS_NUMERICAL_H_

extern float
rssringoccs_Newton_Raphson_Float(float x, float (*f)(float),
                                 float (*f_prime)(float));

extern double
rssringoccs_Newton_Raphson_Double(double x, double (*f)(double),
                                  double (*f_prime)(double));

extern long double
rssringoccs_Newton_Raphson_LongDouble(long double x,
                                      long double (*f)(long double),
                                      long double (*f_prime)(long double));

extern rssringoccs_ComplexDouble
rssringoccs_Newton_Raphson_Complex(
    rssringoccs_ComplexDouble z,
    rssringoccs_ComplexDouble (*f)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_prime)(rssringoccs_ComplexDouble)
);

#endif
/*  End of include guard: #ifndef _RSS_RINGOCCS_NUMERICAL_H_                  */
