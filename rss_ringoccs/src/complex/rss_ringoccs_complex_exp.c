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

/******************************************************************************
 *  Example:                                                                  *
 *      Let's compute the complex exponential of the values 1, i pi, and 1+i. *
 *                                                                            *
 *      #include <rss_ringoccs/src/complex/rss_ringoccs_complex.h>            *
 *      #include <stdio.h>                                                    *
 *                                                                            *
 *      int main(void)                                                        *
 *      {                                                                     *
 *          rssringoccs_ComplexDouble z0, z1, z2;                             *
 *          rssringoccs_ComplexDouble exp_z0, exp_z1, exp_z2;                 *
 *          double re, im, exp_re, exp_im;                                    *
 *                                                                            *
 *          z0 = rssringoccs_Complex_One;                                     *
 *          z1 = rssringoccs_Complex_Rect(0.0, 3.1415926);                    *
 *          z2 = rssringoccs_Complex_Rect(1.0, 1.0);                          *
 *                                                                            *
 *          exp_z0 = rssringoccs_Complex_Exp(z0);                             *
 *          exp_z1 = rssringoccs_Complex_Exp(z1);                             *
 *          exp_z2 = rssringoccs_Complex_Exp(z2);                             *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z0);                           *
 *          im = rssringoccs_Complex_Imag_Part(z0);                           *
 *          exp_re = rssringoccs_Complex_Real_Part(exp_z0);                   *
 *          exp_im = rssringoccs_Complex_Imag_Part(exp_z0);                   *
 *          printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z1);                           *
 *          im = rssringoccs_Complex_Imag_Part(z1);                           *
 *          exp_re = rssringoccs_Complex_Real_Part(exp_z1);                   *
 *          exp_im = rssringoccs_Complex_Imag_Part(exp_z1);                   *
 *          printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z2);                           *
 *          im = rssringoccs_Complex_Imag_Part(z2);                           *
 *          exp_re = rssringoccs_Complex_Real_Part(exp_z2);                   *
 *          exp_im = rssringoccs_Complex_Imag_Part(exp_z2);                   *
 *          printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);     *
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
 *          exp(1.000000 + i0.000000) = 2.718282 + i0.000000                  *
 *          exp(0.000000 + i3.141593) = -1.000000 + i0.000000                 *
 *          exp(1.000000 + i1.000000) = 1.468694 + i2.287355                  *
 *      In agreement with known values of the complex exponential.            *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the complex exponential of a complex number z = x + iy.           */
rssringoccs_ComplexDouble rssringoccs_Complex_Exp(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble exp_z;
    double real, imag;
    double exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_Exp_Double(real);
    exp_z_real = exp_real * rssringoccs_Cos_Double(imag);
    exp_z_imag = exp_real * rssringoccs_Sin_Double(imag);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    exp_z = rssringoccs_Complex_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}
