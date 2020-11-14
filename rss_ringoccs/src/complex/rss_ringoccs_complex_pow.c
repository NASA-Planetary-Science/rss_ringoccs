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
 *                        rss_ringoccs_complex_pow                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex power function. If           *
 *      rss_ringoccs was built with C99 complex.h support, this function is   *
 *      not compiled and instead rssringoccs_Complex_Pow is just an alias for *
 *      the cpow function. By default rss_ringoccs builds with C89 (commonly  *
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
 *          as the prototype for rssringoccs_Complex_Pow.                     *
 ******************************************************************************
 *                                 WARNINGS                                   *
 *  1.) This function implicitly uses the complex log function, and hence     *
 *      there is a branch cut in the second variable. The function is         *
 *      continuous in the first variable, i.e. no branch cut.                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************
 *                                History                                     *
 *  2020/11/14 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/******************************************************************************
 *  Example:                                                                  *
 *      Let's compute the complex power i^i and ((1+i)/sqrt(2))^2             *
 *                                                                            *
 *          #include <rss_ringoccs/src/complex/rss_ringoccs_complex.h>        *
 *          #include <stdio.h>                                                *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              rssringoccs_ComplexDouble z, pow;                             *
 *              double re, im;                                                *
 *                                                                            *
 *              pow = rssringoccs_Complex_Pow(rssringoccs_Imaginary_Unit,     *
 *                                            rssringoccs_Imaginary_Unit);    *
 *              re = rssringoccs_Complex_Real_Part(pow);                      *
 *              im = rssringoccs_Complex_Imag_Part(pow);                      *
 *              printf("i^i = %f + i%f\n", re, im);                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(0.7071067811865476,              *
 *                                           0.7071067811865476);             *
 *              pow = rssringoccs_Complex_Real_Pow(z, 2);                     *
 *              re = rssringoccs_Complex_Real_Part(pow);                      *
 *              im = rssringoccs_Complex_Imag_Part(pow);                      *
 *              printf("((1+i)/sqrt(2))^2 = %f + i%f\n", re, im);             *
 *              return 0;                                                     *
 *          }                                                                 *
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
 *          i^i = 0.207880 + i0.000000                                        *
 *          ((1+i)/sqrt(2))^2 = 0.000000 + i1.000000                          *
 *      So, oddly enough, i^i is a real number. It is exp(-pi/2). This also   *
 *      shows that (1+i)/sqrt(2) is a square root of i.                       *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the complex power z0^z1.                                          */
rssringoccs_ComplexDouble
rssringoccs_Complex_Pow(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z0_to_the_z1, ln_z0, z1_ln_z0;

    /*  We can write x^y as exp(y ln(x)) and this is how we'll compute for    *
     *  complex powers. First compute log(z1).                                */
    ln_z0 = rssringoccs_Complex_Log(z0);

    /*  Next use rssringoccs_Complex_Multiply to compute the product with z0. */
    z1_ln_z0 = rssringoccs_Complex_Multiply(z1, ln_z0);

    /*  And finally exponentiate.                                             */
    z0_to_the_z1 = rssringoccs_Complex_Exp(z1_ln_z0);
    return z0_to_the_z1;
}

/*  Compute the complex power z^x for x real.                                 */
rssringoccs_ComplexDouble
rssringoccs_Complex_Real_Pow(rssringoccs_ComplexDouble z, double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z_to_the_x, ln_z, x_ln_z;

    /*  We can write x^y as exp(y ln(x)) and this is how we'll compute for    *
     *  complex powers. First compute log(z1).                                */
    ln_z = rssringoccs_Complex_Log(z);

    /*  Next use rssringoccs_Complex_Scale to compute the product with x.     */
    x_ln_z = rssringoccs_Complex_Scale(x, ln_z);

    /*  And finally exponentiate.                                             */
    z_to_the_x = rssringoccs_Complex_Exp(x_ln_z);
    return z_to_the_x;
}
