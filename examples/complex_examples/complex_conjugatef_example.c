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
 *      Provides an example of using the complex conjugate function.          *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  Complex functions defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for computing the complex conjugate of various test values.       */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat z[7];
    rssringoccs_ComplexFloat w[7];
    float re_z, im_z, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set the test values in the array z.                                   */
    z[0] = rssringoccs_CFloat_Zero;
    z[1] = rssringoccs_CFloat_One;
    z[2] = rssringoccs_CFloat_I;
    z[3] = rssringoccs_CFloat_Rect(1.0F, 1.0F);
    z[4] = rssringoccs_CFloat_NaN;
    z[5] = rssringoccs_CFloat_Infinity;

    /*  Loop over the results and print them.                                 */
    for (n=0; n<6; ++n)
    {
        /*  Compute the complex conjugate of the nth value.                   */
        w[n] = rssringoccs_CFloat_Conjugate(z[n]);

        /*  Extract the real and imaginary parts from z[n] and w[n].          */
        re_z = rssringoccs_CFloat_Real_Part(z[n]);
        im_z = rssringoccs_CFloat_Imag_Part(z[n]);
        re_w = rssringoccs_CFloat_Real_Part(w[n]);
        im_w = rssringoccs_CFloat_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("conj(%f + i%f) = %f + i%f\n", re_z, im_z, re_w, im_w);
    }
    /*  End of for loop computing conj(z).                                    */

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_conjugatef_example.c -o test -lrssringoccs                *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_conjugatef_example.c -o test -lrssringoccs            *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_conjugatef_example.c -o test -lrssringoccs     *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      conj(0.000000 + i0.000000) = 0.000000 + i-0.000000                    *
 *      conj(1.000000 + i0.000000) = 1.000000 + i-0.000000                    *
 *      conj(0.000000 + i1.000000) = 0.000000 + i-1.000000                    *
 *      conj(1.000000 + i1.000000) = 1.000000 + i-1.000000                    *
 *      conj(nan + inan) = nan + inan                                         *
 *      conj(inf + iinf) = inf + i-inf                                        *
 ******************************************************************************/
