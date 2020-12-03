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
 *                    rss_ringoccs_complex_conjugate                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex conjugate.                   *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************
 *  2020/12/02 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  Single precision complex conjugate function (conjf equivalent).           */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Conjugate(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float real, imag;
    rssringoccs_ComplexFloat conj_z;

    /*  Extract the values from the complex number.                           */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  The complex conjugate of x+iy is just x-iy, compute this.             */
    conj_z = rssringoccs_CFloat_Rect(real, -imag);
    return conj_z;
}
/*  End of rssringoccs_CFloat_Conjugate.                                */

/*  Double precision complex conjugate function (conj equivalent).            */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Conjugate(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag;
    rssringoccs_ComplexDouble conj_z;

    /*  Extract the values from the complex number.                           */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  The complex conjugate of x+iy is just x-iy, compute this.             */
    conj_z = rssringoccs_CDouble_Rect(real, -imag);
    return conj_z;
}
/*  End of rssringoccs_CDouble_Conjugate.                               */

/*  Long double precision complex conjugate function (conjl equivalent).      */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Conjugate(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double real, imag;
    rssringoccs_ComplexLongDouble conj_z;

    /*  Extract the values from the complex number.                           */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  The complex conjugate of x+iy is just x-iy, compute this.             */
    conj_z = rssringoccs_CLDouble_Rect(real, -imag);
    return conj_z;
}
/*  End of rssringoccs_CLDouble_Conjugate.                           */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just alias the          *
 *  functions found in the library.                                           */

/*  Single precision complex conjugate.                                       */
float rssringoccs_CFloat_Conjugate(rssringoccs_ComplexFloat z)
{
    return conjf(z);
}
/*  End of rssringoccs_CFloat_Conjugate.                                */

/*  Double precision complex conjugate.                                       */
double rssringoccs_CDouble_Conjugate(rssringoccs_ComplexDouble z)
{
    return conj(z);
}
/*  End of rssringoccs_CDouble_Conjugate.                               */

/*  Long double precision complex conjugate.                                  */
long double
rssringoccs_CLDouble_Conjugate(rssringoccs_ComplexLongDouble z)
{
    return conjl(z);
}
/*  End of rssringoccs_CLDouble_Conjugate.                           */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
