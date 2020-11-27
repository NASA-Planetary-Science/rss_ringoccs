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
 *                        rss_ringoccs_complex_log                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex log function. If             *
 *      rss_ringoccs was built with C99 complex.h support, this function is   *
 *      not compiled and instead rssringoccs_Complex_Log is just an alias for *
 *      the clog function. By default rss_ringoccs builds with C89 (commonly  *
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
 *          as the prototype for rssringoccs_Complex_Log.                     *
 ******************************************************************************
 *                                 WARNINGS                                   *
 *  1.) This function implicitly uses the atan2 function via the              *
 *      rssringoccs_Complex_Argument function. Because of this the branch cut *
 *      for the complex log occurs along the negative x-axis. No option is    *
 *      provided to choose different branches. One can artificially change    *
 *      the branch by adding a multiple of 2 pi to the imaginary part.        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************/

/******************************************************************************
 *  Example:                                                                  *
 *      Let's compute the complex log of the values 1, -1, and i.             *
 *                                                                            *
 *      #include <rss_ringoccs/src/complex/rss_ringoccs_complex.h>            *
 *      #include <stdio.h>                                                    *
 *                                                                            *
 *      int main(void)                                                        *
 *      {                                                                     *
 *          rssringoccs_ComplexDouble z0, z1, z2;                             *
 *          rssringoccs_ComplexDouble log_z0, log_z1, log_z2;                 *
 *          double re, im, log_re, log_im;                                    *
 *                                                                            *
 *          z0 = rssringoccs_Complex_One;                                     *
 *          z1 = rssringoccs_Complex_Rect(-1.0, 0.0);                         *
 *          z2 = rssringoccs_Imaginary_Unit;                                  *
 *                                                                            *
 *          log_z0 = rssringoccs_Complex_Log(z0);                             *
 *          log_z1 = rssringoccs_Complex_Log(z1);                             *
 *          log_z2 = rssringoccs_Complex_Log(z2);                             *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z0);                           *
 *          im = rssringoccs_Complex_Imag_Part(z0);                           *
 *          log_re = rssringoccs_Complex_Real_Part(log_z0);                   *
 *          log_im = rssringoccs_Complex_Imag_Part(log_z0);                   *
 *          printf("log(%f + i%f) = %f + i%f\n", re, im, log_re, log_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z1);                           *
 *          im = rssringoccs_Complex_Imag_Part(z1);                           *
 *          log_re = rssringoccs_Complex_Real_Part(log_z1);                   *
 *          log_im = rssringoccs_Complex_Imag_Part(log_z1);                   *
 *          printf("log(%f + i%f) = %f + i%f\n", re, im, log_re, log_im);     *
 *                                                                            *
 *          re = rssringoccs_Complex_Real_Part(z2);                           *
 *          im = rssringoccs_Complex_Imag_Part(z2);                           *
 *          log_re = rssringoccs_Complex_Real_Part(log_z2);                   *
 *          log_im = rssringoccs_Complex_Imag_Part(log_z2);                   *
 *          printf("log(%f + i%f) = %f + i%f\n", re, im, log_re, log_im);     *
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
 *          log(1.000000 + i0.000000) = 0.000000 + i0.000000                  *
 *          log(-1.000000 + i0.000000) = 0.000000 + i3.141593                 *
 *          log(0.000000 + i1.000000) = 0.000000 + i1.570796                  *
 *      In agreement with known values of the complex logarithm. That is,     *
 *      log(1) = 0, log(-1) = ipi, and log(i) = i pi / 2.                     *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the complex log of a complex number z = r exp(i theta) where      *
 *  theta is a real number between -pi and pi.                                */
rssringoccs_ComplexDouble
rssringoccs_Complex_Log(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double r, theta, real;
    rssringoccs_ComplexDouble ln_z;

    /*  Get the polar representation of the complex number z.                 */
    r = rssringoccs_Complex_Abs(z);
    theta = rssringoccs_Complex_Argument(z);

    /*  The real part is just ln(r), and the imaginary part is theta.         */
    real = rssringoccs_Double_Log(r);

    /*  Use rssringoccs_Complex_Rect to create the complex number and return. */
    ln_z = rssringoccs_Complex_Rect(real, theta);
    return ln_z;
}
