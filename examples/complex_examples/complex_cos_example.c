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
 ******************************************************************************/

/*  Let's compute the complex cosine of the values pi, i pi, and 0.           */

/*  If rss_ringoccs built correctly, rss_ringoccs_complex.h is located in     *
 *  /usr/local/include/rss_ringoccs/include. We can include this as follows:  */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  ONE_PI is defined here.                                                   */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

int main(void)
{
    /*  Declare all necessary variables. C89 requires declarations at the top *
     *  of a block.                                                           */
    rssringoccs_ComplexDouble z[3];
    rssringoccs_ComplexDouble w[3];
    double re_z, im_z, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set z0, z1, and z2 to 0, i pi, and pi, respectively.                  */
    z[0] = rssringoccs_Complex_Zero;
    z[1] = rssringoccs_Complex_Rect(0.0, ONE_PI);
    z[2] = rssringoccs_Complex_Rect(ONE_PI, 0.0);

    /*  Loop over the results and print them.                                 */
    for (n=0; n<3; ++n)
    {
        /*  Compute the complex cosine of the nth value.                      */
        w[n] = rssringoccs_Complex_Cos(z[n]);

        /*  Extract the real and imaginary parts from z[n] and w[n].          */
        re_z = rssringoccs_Complex_Real_Part(z[n]);
        im_z = rssringoccs_Complex_Imag_Part(z[n]);
        re_w = rssringoccs_Complex_Real_Part(w[n]);
        im_w = rssringoccs_Complex_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("cos(%f + i%f) = %f + i%f\n", re_z, im_z, re_w, im_w);
    }

    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_cos_example.c -o test -lrssringoccs                       *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_src.sh), then change       *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_cos_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      cos(0.000000 + i0.000000) = 1.000000 + i-0.000000                     *
 *      cos(0.000000 + i3.141593) = 11.591953 + i-0.000000                    *
 *      cos(3.141593 + i0.000000) = -1.000000 + i-0.000000                    *
 *  In agreement with known values of the complex cosine.                     *
 ******************************************************************************/
