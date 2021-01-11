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
 *                        rss_ringoccs_complex_add                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex addition.                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  This function is equivalent to the creal function in complex.h (C99).     */
float rssringoccs_CFloat_Real_Part(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float real;

    /*  The real component is stored as the first entry in the dat array      *
     *  contained in a rssringoccs_ComplexDouble struct. Return this.         */
    real = z.dat[0];
    return real;
}

double rssringoccs_CDouble_Real_Part(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real;

    /*  The real component is stored as the first entry in the dat array      *
     *  contained in a rssringoccs_ComplexDouble struct. Return this.         */
    real = z.dat[0];
    return real;
}

long double
rssringoccs_CLDouble_Real_Part(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double real;

    /*  The real component is stored as the first entry in the dat array      *
     *  contained in a rssringoccs_ComplexDouble struct. Return this.         */
    real = z.dat[0];
    return real;
}

#else

double rssringoccs_CFloat_Real_Part(rssringoccs_ComplexFloat z)
{
    return crealf(z);
}

double rssringoccs_CDouble_Real_Part(rssringoccs_ComplexDouble z)
{
    return creal(z);
}

long double
rssringoccs_CLDouble_Real_Part(rssringoccs_ComplexLongDouble z)
{
    return creall(z);
}

#endif
