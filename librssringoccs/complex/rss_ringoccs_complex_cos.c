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
 *                        rss_ringoccs_complex_cos                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex cosine function.             *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header file where rssringoccs_ComplexDouble is defined, as well   *
 *          as the prototype for rssringoccs_Complex_Cos.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/11/14 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  Single precision complex cosine.                                          */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Cos(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float x, y, real, imag;
    rssringoccs_ComplexFloat cos_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CFloat_Real_Part(z);
    y = rssringoccs_CFloat_Imag_Part(z);

    /*  The real part is cos(x)cosh(y).                                       */
    real = rssringoccs_Float_Cos(x) * rssringoccs_Float_Cosh(y);

    /*  And the imaginary part is -sin(x)sinh(y).                             */
    imag = -rssringoccs_Float_Sin(x) * rssringoccs_Float_Sinh(y);

    /*  Create the output and return.                                         */
    cos_z = rssringoccs_CFloat_Rect(real, imag);
    return cos_z;
}
/*  End of rssringoccs_CFloat_Cos.                                            */

/*  Double precision complex cosine.                                          */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Cos(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, real, imag;
    rssringoccs_ComplexDouble cos_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CDouble_Real_Part(z);
    y = rssringoccs_CDouble_Imag_Part(z);

    /*  The real part is cos(x)cosh(y).                                       */
    real = rssringoccs_Double_Cos(x) * rssringoccs_Double_Cosh(y);

    /*  And the imaginary part is -sin(x)sinh(y).                             */
    imag = -rssringoccs_Double_Sin(x) * rssringoccs_Double_Sinh(y);

    /*  Create the output and return.                                         */
    cos_z = rssringoccs_CDouble_Rect(real, imag);
    return cos_z;
}
/*  End of rssringoccs_CDouble_Cos.                                           */

/*  Long double precision complex cosine.                                     */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Cos(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double x, y, real, imag;
    rssringoccs_ComplexLongDouble cos_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CLDouble_Real_Part(z);
    y = rssringoccs_CLDouble_Imag_Part(z);

    /*  The real part is cos(x)cosh(y).                                       */
    real = rssringoccs_LDouble_Cos(x) * rssringoccs_LDouble_Cosh(y);

    /*  And the imaginary part is -sin(x)sinh(y).                             */
    imag = -rssringoccs_LDouble_Sin(x) * rssringoccs_LDouble_Sinh(y);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    cos_z = rssringoccs_CLDouble_Rect(real, imag);
    return cos_z;
}
/*  End of rssringoccs_CLDouble_Cos.                                 */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just alias the          *
 *  functions found in the library.                                           */

/*  Single precision complex cosine.                                          */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Cos(rssringoccs_ComplexFloat z)
{
    return ccosf(z);
}
/*  End of rssringoccs_CLDouble_Cos.                                 */

/*  Double precision complex cosine.                                          */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Cos(rssringoccs_ComplexDouble z)
{
    return ccos(z);
}
/*  End of rssringoccs_CDouble_Cos.                                     */

/*  Long double precision complex cosine.                                     */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Cos(rssringoccs_ComplexLongDouble z)
{
    return ccosl(z);
}
/*  End of rssringoccs_CLDouble_Cos.                                 */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
