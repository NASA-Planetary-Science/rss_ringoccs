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
 *      Provides an example of using the complex exponential function.        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  Let's compute the complex exponential of the values 1, i, i pi, and pi.   */

/*  Complex exponential is declared here.                                     */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  rssringoccs_One_Pi is defined here.                                       */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for computing the complex exponential of 1, i, i pi, and pi.      */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z[4];
    rssringoccs_ComplexDouble w[4];
    double re_z, im_z, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set z0, z1, z2, and z3 to 1, i, i pi, and pi, respectively.           */
    z[0] = rssringoccs_CDouble_One;
    z[1] = rssringoccs_CDouble_I;
    z[2] = rssringoccs_CDouble_Rect(0.0, rssringoccs_One_Pi);
    z[3] = rssringoccs_CDouble_Rect(rssringoccs_One_Pi, 0.0);

    /*  Loop over the results and print them.                                 */
    for (n=0; n<4; ++n)
    {
        /*  Compute the complex exponential of the nth value.                 */
        w[n] = rssringoccs_CDouble_Exp(z[n]);

        /*  Extract the real and imaginary parts from z[n] and w[n].          */
        re_z = rssringoccs_CDouble_Real_Part(z[n]);
        im_z = rssringoccs_CDouble_Imag_Part(z[n]);
        re_w = rssringoccs_CDouble_Real_Part(w[n]);
        im_w = rssringoccs_CDouble_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("exp(%f + i%f) = %f + i%f\n", re_z, im_z, re_w, im_w);
    }

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_exp_example.c -o test -lrssringoccs                       *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_exp_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_exp_example.c -o test -lrssringoccs            *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      exp(1.000000 + i0.000000) = 2.718282 + i0.000000                      *
 *      exp(0.000000 + i1.000000) = 0.540302 + i0.841471                      *
 *      exp(0.000000 + i3.141593) = -1.000000 + i0.000000                     *
 *      exp(3.141593 + i0.000000) = 23.140693 + i0.000000                     *
 *  In agreement with known values of the complex exponential.                *
 ******************************************************************************/
