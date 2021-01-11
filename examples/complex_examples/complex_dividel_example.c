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
 *      Provides an example of dividing complex numbers.                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 22, 2020                                             *
 ******************************************************************************/

/*  Complex functions defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Real nan and inf found here.                                              */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for dividing two complex numbers.                                 */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble z0[7], z1[7], w[7];
    long double re_z0, im_z0, re_z1, im_z1, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set the test values in the array z0.                                  */
    z0[0] = rssringoccs_CLDouble_Zero;
    z0[1] = rssringoccs_CLDouble_One;
    z0[2] = rssringoccs_CLDouble_Rect(1.0L, 1.0L);
    z0[3] = rssringoccs_CLDouble_Rect(rssringoccs_NaN_L, 0.0L);
    z0[4] = rssringoccs_CLDouble_Rect(rssringoccs_Infinity_L, 0.0L);
    z0[5] = rssringoccs_CLDouble_NaN;
    z0[6] = rssringoccs_CLDouble_Infinity;

    /*  Set the test values for the array z1.                                 */
    z1[0] = rssringoccs_CLDouble_Rect(1.0L, 1.0L);
    z1[1] = rssringoccs_CLDouble_Rect(2.0L, 0.0L);
    z1[2] = rssringoccs_CLDouble_Rect(1.0L, 3.0L);
    z1[3] = rssringoccs_CLDouble_Rect(1.0L, 4.0L);
    z1[4] = rssringoccs_CLDouble_Rect(1.0L, 0.0L);
    z1[5] = rssringoccs_CLDouble_Rect(1.0L, 8.0L);
    z1[6] = rssringoccs_CLDouble_Rect(1.0L, 1.0L);

    /*  Loop over the results and print them.                                 */
    for (n=0; n<7; ++n)
    {
        /*  Compute z0 / z1 of the nth value.                                 */
        w[n] = rssringoccs_CLDouble_Divide(z0[n], z1[n]);

        /*  Extract the real and imaginary parts from z0[n].                  */
        re_z0 = rssringoccs_CLDouble_Real_Part(z0[n]);
        im_z0 = rssringoccs_CLDouble_Imag_Part(z0[n]);

        /*  Extract the real and imaginary parts from z1[n].                  */
        re_z1 = rssringoccs_CLDouble_Real_Part(z1[n]);
        im_z1 = rssringoccs_CLDouble_Imag_Part(z1[n]);

        /*  Extract the real and imaginary parts from w[n].                   */
        re_w = rssringoccs_CLDouble_Real_Part(w[n]);
        im_w = rssringoccs_CLDouble_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("(%Lf + i%Lf) / (%Lf + i%Lf) = %Lf +i%Lf\n",
               re_z0, im_z0, re_z1, im_z1, re_w, im_w);
    }
    /*  End of for loop z0 / z1.                                              */

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_dividel_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_dividel_example.c -o test -lrssringoccs               *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_dividel_example.c -o test -lrssringoccs        *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      (0.000000 + i0.000000) / (1.000000 + i1.000000) = 0.000000 +i0.000000 *
 *      (1.000000 + i0.000000) / (2.000000 + i0.000000) = 0.500000 +i0.000000 *
 *      (1.000000 + i1.000000) / (1.000000 + i3.000000) = 0.400000 +i-0.200000*
 *      (nan + i0.000000) / (1.000000 + i4.000000) = nan +inan                *
 *      (inf + i0.000000) / (1.000000 + i0.000000) = inf +inan                *
 *      (nan + inan) / (1.000000 + i8.000000) = nan +inan                     *
 *      (inf + iinf) / (1.000000 + i1.000000) = inf +inan                     *
 ******************************************************************************/
