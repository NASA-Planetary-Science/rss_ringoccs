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
 *                        rss_ringoccs_complex_abs                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex modulus (absolute value).        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/01 (Ryan Maguire):                                                *
 *      Added abs squared functions.                                          *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file which contains aliases for the functions in the standard C    *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  Single precision complex abs function (cabsf equivalent).                 */
float rssringoccs_ComplexFloat_Abs(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float real, imag, abs_value;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexFloat_Real_Part(z);
    imag = rssringoccs_ComplexFloat_Imag_Part(z);

    /*  The absolute value is just sqrt(x^2 + y^2) so compute this.           */
    abs_value = rssringoccs_Float_Sqrt(real*real + imag*imag);
    return abs_value;
}
/*  End of rssringoccs_ComplexFloat_Abs.                                      */

/*  Double precision complex abs function (cabs equivalent).                  */
double rssringoccs_ComplexDouble_Abs(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag, abs_value;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexDouble_Real_Part(z);
    imag = rssringoccs_ComplexDouble_Imag_Part(z);

    /*  The absolute value is just sqrt(x^2 + y^2) so compute this.           */
    abs_value = rssringoccs_Double_Sqrt(real*real + imag*imag);
    return abs_value;
}
/*  End of rssringoccs_ComplexDouble_Abs.                                     */

/*  Long double precision complex abs function (cabsl equivalent).            */
long double rssringoccs_ComplexLongDouble_Abs(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double real, imag, abs_value;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexLongDouble_Real_Part(z);
    imag = rssringoccs_ComplexLongDouble_Imag_Part(z);

    /*  The absolute value is just sqrt(x^2 + y^2) so compute this.           */
    abs_value = rssringoccs_LongDouble_Sqrt(real*real + imag*imag);
    return abs_value;
}
/*  End of rssringoccs_ComplexLongDouble_Abs.                                 */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just alias the          *
 *  functions found in the library.                                           */

/*  Single precision absolute value function, alias for cabsf.                */
float rssringoccs_ComplexFloat_Abs(rssringoccs_ComplexFloat z)
{
    return cabsf(z);
}
/*  End of rssringoccs_ComplexFloat_Abs.                                      */

/*  Double precision absolute value function, alias for cabs.                 */
double rssringoccs_ComplexDouble_Abs(rssringoccs_ComplexDouble z)
{
    return cabs(z);
}
/*  End of rssringoccs_ComplexDouble_Abs.                                     */

/*  Long double precision absolute value function, alias for cabsl.           */
long double rssringoccs_ComplexLongDouble_Abs(rssringoccs_ComplexLongDouble z)
{
    return cabsl(z);
}
/*  End of rssringoccs_ComplexLongDouble_Abs.                                 */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */

/*  C99 does not provide |z|^2 functions, so we can use the same algorithm    *
 *  for both the C89 and C99 renditions of rss_ringoccs.                      */

/*  Single precision abs squared function.                                    */
float rssringoccs_ComplexFloat_Abs_Squared(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float real, imag, abs_sq;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexFloat_Real_Part(z);
    imag = rssringoccs_ComplexFloat_Imag_Part(z);

    /*  |z|^2 = x^2 + y^2 so compute this.                                    */
    abs_sq = real*real + imag*imag;
    return abs_sq;
}
/*  End of rssringoccs_ComplexFloat_Abs_Square.                               */

/*  Double precision abs squared function.                                    */
double rssringoccs_ComplexDouble_Abs_Squared(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag, abs_sq;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexDouble_Real_Part(z);
    imag = rssringoccs_ComplexDouble_Imag_Part(z);

    /*  |z|^2 = x^2 + y^2 so compute this.                                    */
    abs_sq = real*real + imag*imag;
    return abs_sq;
}
/*  End of rssringoccs_ComplexDouble_Abs_Square.                              */

/*  Long double precision abs squared function.                               */
long double
rssringoccs_ComplexLongDouble_Abs_Squared(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double real, imag, abs_sq;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_ComplexLongDouble_Real_Part(z);
    imag = rssringoccs_ComplexLongDouble_Imag_Part(z);

    /*  |z|^2 = x^2 + y^2 so compute this.                                    */
    abs_sq = real*real + imag*imag;
    return abs_sq;
}
/*  End of rssringoccs_ComplexLongDouble_Abs_Square.                          */
