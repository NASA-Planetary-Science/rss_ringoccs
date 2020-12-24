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
 *      Provides an example of using the complex cosine function.             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 23, 2020                                             *
 ******************************************************************************/

/*  Complex cosine is declared here.                                          */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  rssringoccs_One_Pi and real nan and inf are defined here.                 */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for computing the complex cosine of various test values.          */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble z[7];
    rssringoccs_ComplexLongDouble w[7];
    long double re_z, im_z, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set the test values in the array z.                                   */
    z[0] = rssringoccs_CLDouble_Zero;
    z[1] = rssringoccs_CLDouble_Rect(0.0L, rssringoccs_One_Pi_L);
    z[2] = rssringoccs_CLDouble_Rect(rssringoccs_One_Pi_L, 0.0L);
    z[3] = rssringoccs_CLDouble_Rect(rssringoccs_NaN_L, 0.0L);
    z[4] = rssringoccs_CLDouble_Rect(rssringoccs_Infinity_L, 0.0L);
    z[5] = rssringoccs_CLDouble_NaN;
    z[6] = rssringoccs_CLDouble_Infinity;

    /*  Loop over the results and print them.                                 */
    for (n=0; n<7; ++n)
    {
        /*  Compute the complex cosine of the nth value.                      */
        w[n] = rssringoccs_CLDouble_Cos(z[n]);

        /*  Extract the real and imaginary parts from z[n] and w[n].          */
        re_z = rssringoccs_CLDouble_Real_Part(z[n]);
        im_z = rssringoccs_CLDouble_Imag_Part(z[n]);
        re_w = rssringoccs_CLDouble_Real_Part(w[n]);
        im_w = rssringoccs_CLDouble_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("cos(%Lf + i%Lf) = %Lf + i%Lf\n", re_z, im_z, re_w, im_w);
    }
    /*  End of for loop computing cos(z).                                     */

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_cosl_example.c -o test -lrssringoccs                      *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_cosl_example.c -o test -lrssringoccs                  *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_cosl_example.c -o test -lrssringoccs           *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      cos(0.000000 + i0.000000) = 1.000000 + i-0.000000                     *
 *      cos(0.000000 + i3.141593) = 11.591953 + i-0.000000                    *
 *      cos(3.141593 + i0.000000) = -1.000000 + i-0.000000                    *
 *      cos(nan + i0.000000) = nan + inan                                     *
 *      cos(inf + i0.000000) = nan + inan                                     *
 *      cos(nan + inan) = nan + inan                                          *
 *      cos(inf + iinf) = nan + inan                                          *
 ******************************************************************************/
