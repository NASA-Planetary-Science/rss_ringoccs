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
 *                        rss_ringoccs_complex_exp                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex exponention function. If     *
 *      rss_ringoccs was built with C99 complex.h support, this function is   *
 *      not compiled and instead rssringoccs_Complex_Exp is just an alias for *
 *      the cexp function. By default rss_ringoccs builds with C89 (commonly  *
 *      called ANSI C) support, so this file will be a part of the build.     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header file where rssringoccs_ComplexDouble is defined, as well   *
 *          as the prototype for rssringoccs_Complex_Exp.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the complex exponential of a complex number z = x + iy.           */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Exp(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble exp_z;
    double real, imag;
    double exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_Double_Exp(real);

    /*  In the case that z is real, use the real valued exponential. This     *
     *  avoid the result of exp(inf) = inf + i nan. The imaginary part of     *
     *  complex exp(inf) will be exp(inf) * sin(0) = inf * 0 which results in *
     *  nan. This if-then statement avoids this.                              */
    if (imag == 0.0)
    {
        exp_z_real = exp_real;
        exp_z_imag = 0.0;
    }

    /*  When we have non-zero imaginary part, resort to Euler's formula.      */
    else
    {
        exp_z_real = exp_real * rssringoccs_Double_Cos(imag);
        exp_z_imag = exp_real * rssringoccs_Double_Sin(imag);
    }

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    exp_z = rssringoccs_CDouble_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}
