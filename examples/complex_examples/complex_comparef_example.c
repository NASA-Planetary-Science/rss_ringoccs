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
 *      Provides an example of using the compare function for complex numbers.*
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 22, 2020                                             *
 ******************************************************************************/

/*  Complex functions defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for testing the compare function.                                 */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat z0[9], z1[9];
    float re_z0, im_z0, re_z1, im_z1;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set the test values in the array z0.                                  */
    z0[0] = rssringoccs_CFloat_Zero;
    z0[1] = rssringoccs_CFloat_One;
    z0[2] = rssringoccs_CFloat_Rect(1.0F, 1.0F);
    z0[3] = rssringoccs_CFloat_I;
    z0[4] = rssringoccs_CFloat_Rect(-1.0F, 1.0F);
    z0[5] = rssringoccs_CFloat_Rect(-1.0F, 0.0F);
    z0[6] = rssringoccs_CFloat_Rect(-1.0F, -1.0F);
    z0[7] = rssringoccs_CFloat_Rect(0.0F, -1.0F);
    z0[8] = rssringoccs_CFloat_Rect(1.0F, -1.0F);

    /*  Set the test values in the array z1.                                  */
    z1[0] = rssringoccs_CFloat_Zero;
    z1[1] = rssringoccs_CFloat_Rect(1.0F, 1.0F);
    z1[2] = rssringoccs_CFloat_One;
    z1[3] = rssringoccs_CFloat_I;
    z1[4] = rssringoccs_CFloat_Rect(-1.0F, 0.0F);
    z1[5] = rssringoccs_CFloat_Rect(-1.0F, 1.0F);
    z1[6] = rssringoccs_CFloat_Rect(-1.0F, -1.0F);
    z1[7] = rssringoccs_CFloat_Rect(1.0F, -1.0F);
    z1[8] = rssringoccs_CFloat_Rect(0.0F, -1.0F);

    /*  Loop over the results and print them.                                 */
    for (n=0; n<9; ++n)
    {
        /*  Extract the real and imaginary parts from z0[n].                  */
        re_z0 = rssringoccs_CFloat_Real_Part(z0[n]);
        im_z0 = rssringoccs_CFloat_Imag_Part(z0[n]);

        /*  Extract the real and imaginary parts from z1[n].                  */
        re_z1 = rssringoccs_CFloat_Real_Part(z1[n]);
        im_z1 = rssringoccs_CFloat_Imag_Part(z1[n]);

        /*  Check if the values are equal.                                    */
        if (rssringoccs_CFloat_Compare(z0[n], z1[n]))
            printf("%f + i%f = %f + i%f: True\n", re_z0, im_z0, re_z1, im_z1);
        else
            printf("%f + i%f = %f + i%f: False\n", re_z0, im_z0, re_z1, im_z1);
    }
    /*  End of for loop comparing z0 and z1.                                  */

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_comparef_example.c -o test -lrssringoccs                  *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_comparef_example.c -o test -lrssringoccs              *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic complex_comparef_example.c -o test -lrssringoccs       *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      0.000000 + i0.000000 = 0.000000 + i0.000000: True                     *
 *      1.000000 + i0.000000 = 1.000000 + i1.000000: False                    *
 *      1.000000 + i1.000000 = 1.000000 + i0.000000: False                    *
 *      0.000000 + i1.000000 = 0.000000 + i1.000000: True                     *
 *      -1.000000 + i1.000000 = -1.000000 + i0.000000: False                  *
 *      -1.000000 + i0.000000 = -1.000000 + i1.000000: False                  *
 *      -1.000000 + i-1.000000 = -1.000000 + i-1.000000: True                 *
 *      0.000000 + i-1.000000 = 1.000000 + i-1.000000: False                  *
 *      1.000000 + i-1.000000 = 0.000000 + i-1.000000: False                  *
 ******************************************************************************/
