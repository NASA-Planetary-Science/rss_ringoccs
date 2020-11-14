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
 *      Contains the source code for the complex cosine function. If          *
 *      rss_ringoccs was built with C99 complex.h support, this function is   *
 *      not compiled and instead rssringoccs_Complex_Cos is just an alias for *
 *      the ccos function. By default rss_ringoccs builds with C89 (commonly  *
 *      called ANSI C) support, so this file will be a part of the build.     *
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

/******************************************************************************
 *  Example:                                                                  *
 *      Let's compute the complex cosine of the values pi, i pi, and 0.       *
 *                                                                            *
 *      #include <rss_ringoccs/rss_ringoccs_complex.h>                        *
 *      #include <stdio.h>                                                    *
 *                                                                            *
 *      int main(void)                                                        *
 *      {                                                                     *
 *          rssringoccs_ComplexDouble z0, z1, z2;                             *
 *          rssringoccs_ComplexDouble cos_z0, cos_z1, cos_z2;                 *
 *          double re, im, cos_re, cos_im;                                    *
 *                                                                            *
 *          z0 = rssringoccs_Complex_Zero;                                    *
 *          z1 = rssringoccs_Complex_Rect(0.0, 3.1415926);                    *
 *          z2 = rssringoccs_Complex_Rect(3.1415926, 0.0);                    *
 *                                                                            *
 *          cos_z0 = rssringoccs_Complex_Cos(z0);                             *
 *          cos_z1 = rssringoccs_Complex_Cos(z1);                             *
 *          cos_z2 = rssringoccs_Complex_Cos(z2);                             *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z0);                           *
 *          im = rssringoccs_Complex_Imag_Part(z0);                           *
 *          cos_re = rssringoccs_Complex_Real_Part(cos_z0);                   *
 *          cos_im = rssringoccs_Complex_Imag_Part(cos_z0);                   *
 *          printf("cos(%f + i%f) = %f + i%f\n", re, im, cos_re, cos_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z1);                           *
 *          im = rssringoccs_Complex_Imag_Part(z1);                           *
 *          cos_re = rssringoccs_Complex_Real_Part(cos_z1);                   *
 *          cos_im = rssringoccs_Complex_Imag_Part(cos_z1);                   *
 *          printf("cos(%f + i%f) = %f + i%f\n", re, im, cos_re, cos_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z2);                           *
 *          im = rssringoccs_Complex_Imag_Part(z2);                           *
 *          cos_re = rssringoccs_Complex_Real_Part(cos_z2);                   *
 *          cos_im = rssringoccs_Complex_Imag_Part(cos_z2);                   *
 *          printf("cos(%f + i%f) = %f + i%f\n", re, im, cos_re, cos_im);     *
 *                                                                            *
 *          return 0;                                                         *
 *      }                                                                     *
 *                                                                            *
 *      Naming this test.c and placing it in rss_ringoccs/src/ we can compile *
 *      this with:                                                            *
 *                                                                            *
 *          gcc -I../../ -L/usr/local/lib/ test.c -o test -lrssringoccs       *
 *                                                                            *
 *      If librssringoccs is not in /usr/local/lib/ (this is the default      *
 *      location it is placed in when built via config_src.sh), then change   *
 *      the -L option to the correct location. Change the -I location so that *
 *      rss_ringoccs/ is in your path, if needed.                             *
 *                                                                            *
 *      Running the executable with ./test, this outputs:                     *
 *          cos(0.000000 + i0.000000) = 1.000000 + i-0.000000                 *
 *          cos(0.000000 + i3.141593) = 11.591953 + i-0.000000                *
 *          cos(3.141593 + i0.000000) = -1.000000 + i-0.000000                *
 *      In agreement with known values of the complex cosine.                 *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the cosine of a complex number.                                   */
rssringoccs_ComplexDouble rssringoccs_Complex_Cos(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, real, imag;
    rssringoccs_ComplexDouble cos_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_Complex_Real_Part(z);
    y = rssringoccs_Complex_Imag_Part(z);

    /*  The real part is cos(x)cosh(y).                                       */
    real = rssringoccs_Cos_Double(x)*rssringoccs_Cosh_Double(y);

    /*  And the imaginary part is -sin(x)sinh(y).                             */
    imag = -rssringoccs_Sin_Double(x)*rssringoccs_Sinh_Double(y);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    cos_z = rssringoccs_Complex_Rect(real, imag);
    return cos_z;
}
