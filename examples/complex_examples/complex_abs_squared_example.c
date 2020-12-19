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
 *      Provides an example of using the abs squared function.                *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 17, 2020                                             *
 ******************************************************************************/

/*  Let's compute |z|^2 for 0, 1, 1+i, nan, infinity, complex nan, and        *
 *  complex infinity.                                                         */

/*  Complex functions defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Real nan and inf found here.                                              */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for computing |z|^2 for a few test values.                        */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z[7];
    double w[7];
    double re_z, im_z;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set the test values in the array z.                                   */
    z[0] = rssringoccs_CDouble_Zero;
    z[1] = rssringoccs_CDouble_One;
    z[2] = rssringoccs_CDouble_Rect(1.0, 1.0);
    z[3] = rssringoccs_CDouble_Rect(rssringoccs_NaN, 0.0);
    z[4] = rssringoccs_CDouble_Rect(rssringoccs_Infinity, 0.0);
    z[5] = rssringoccs_CDouble_NaN;
    z[6] = rssringoccs_CDouble_Infinity;

    /*  Loop over the results and print them.                                 */
    for (n=0; n<7; ++n)
    {
        /*  Compute |z|^2 of the nth value.                                   */
        w[n] = rssringoccs_CDouble_Abs_Squared(z[n]);

        /*  Extract the real and imaginary parts from z[n].                   */
        re_z = rssringoccs_CDouble_Real_Part(z[n]);
        im_z = rssringoccs_CDouble_Imag_Part(z[n]);

        /*  And finally, print the result to the screen.                      */
        printf("|%f + i%f|^2 = %f\n", re_z, im_z, w[n]);
    }
    /*  End of for loop computing |z|^2.                                      */

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_abs_squared_example.c -o test -lrssringoccs               *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_abs_squared_example.c -o test -lrssringoccs           *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_abs_squared_example.c -o test -lrssringoccs    *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      |0.000000 + i0.000000|^2 = 0.000000                                   *
 *      |1.000000 + i0.000000|^2 = 1.000000                                   *
 *      |1.000000 + i1.000000|^2 = 2.000000                                   *
 *      |nan + i0.000000|^2 = nan                                             *
 *      |inf + i0.000000|^2 = inf                                             *
 *      |nan + inan|^2 = nan                                                  *
 *      |inf + iinf|^2 = inf                                                  *
 ******************************************************************************/
