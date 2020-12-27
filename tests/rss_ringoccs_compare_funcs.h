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
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide functions for comparing the accuracy and efficiency of        *
 *      functions in rss_ringoccs as opposed to other libraries.              *
 *  NOTE:                                                                     *
 *      librssringoccs does not have any dependencies and will compile on any *
 *      compiler capable of handling C89/C90 or C99 compliant code. The tests *
 *      using these functions use external libraries to compare the results   *
 *      of rss_ringoccs with others. To run these tests requires having these *
 *      libraries available. These tests are NOT required to use rss_ringoccs *
 *      and are mainly for internal use.                                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 12, 2020                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_COMPARE_FUNCS_H__
#define __RSS_RINGOCCS_COMPARE_FUNCS_H__

#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <complex.h>

extern void
rssringoccs_Compare_Float_Funcs(const char *f0_name, float (*f0)(float),
                                const char *f1_name, float (*f1)(float),
                                float start, float end, unsigned long N);

extern void
rssringoccs_Compare_Double_Funcs(const char *f0_name, double (*f0)(double),
                                 const char *f1_name, double (*f1)(double),
                                 double start, double end, unsigned long N);


extern void
rssringoccs_Compare_LDouble_Funcs(const char *f0_name,
                                  long double (*f0)(long double),
                                  const char *f1_name,
                                  long double (*f1)(long double),
                                  long double start, long double end,
                                  unsigned long N);

extern void
rssringoccs_Compare_CFloat_Funcs(const char *f0_name,
                                 rssringoccs_ComplexFloat
                                   (*f0)(rssringoccs_ComplexFloat),
                                 const char *f1_name,
                                 float _Complex
                                   (*f1)(float _Complex),
                                 float start, float end, unsigned long N);

extern void
rssringoccs_Compare_CDouble_Funcs(
    const char *f0_name,
    rssringoccs_ComplexDouble (*f0)(rssringoccs_ComplexDouble),
    const char *f1_name,
    double _Complex (*f1)(double _Complex),
    const double start,
    const double end,
    const unsigned long N);

extern void
rssringoccs_Compare_CLDouble_Funcs(const char *f0_name,
                                   rssringoccs_ComplexLongDouble
                                     (*f0)(rssringoccs_ComplexLongDouble),
                                   const char *f1_name,
                                   long double _Complex
                                     (*f1)(long double _Complex),
                                   long double start, long double end,
                                   unsigned long N);

extern void
rssringoccs_RelCompare_CFloat_Funcs(const char *f0_name,
                                    rssringoccs_ComplexFloat
                                      (*f0)(rssringoccs_ComplexFloat),
                                    const char *f1_name,
                                    float _Complex
                                      (*f1)(float _Complex),
                                    float start, float end, unsigned long N);

extern void
rssringoccs_RelCompare_CDouble_Funcs(const char *f0_name,
                                     rssringoccs_ComplexDouble
                                       (*f0)(rssringoccs_ComplexDouble),
                                     const char *f1_name,
                                     double _Complex
                                       (*f1)(double _Complex),
                                     double start, double end, unsigned long N);

extern void
rssringoccs_RelCompare_CLDouble_Funcs(const char *f0_name,
                                      rssringoccs_ComplexLongDouble
                                        (*f0)(rssringoccs_ComplexLongDouble),
                                      const char *f1_name,
                                      long double _Complex
                                        (*f1)(long double _Complex),
                                      long double start, long double end,
                                      unsigned long N);

extern void
rssringoccs_Compare_Real_CFloat_Funcs(const char *f0_name,
                                      float (*f0)(rssringoccs_ComplexFloat),
                                      const char *f1_name,
                                      float (*f1)(float _Complex),
                                      float start, float end, unsigned long N);

extern void
rssringoccs_Compare_Real_CDouble_Funcs(const char *f0_name,
                                       double (*f0)(rssringoccs_ComplexDouble),
                                       const char *f1_name,
                                       double (*f1)(double _Complex),
                                       double start, double end,
                                       unsigned long N);


extern void
rssringoccs_Compare_Real_CLDouble_Funcs(const char *f0_name,
                                        long double
                                          (*f0)(rssringoccs_ComplexLongDouble),
                                        const char *f1_name,
                                        long double (*f1)(long double _Complex),
                                        long double start, long double end,
                                        unsigned long N);

#endif
/*  End of include guard.                                                     */
