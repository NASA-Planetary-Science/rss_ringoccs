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
 *      Provides an example of using the Euclidean norm function for a three  *
 *      dimensional vector.                                                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  Let's compute the norm of the vector (1, 2, 3).                           */

/*  The Euclidean norm is declared here.                                      */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Routine for computing the norm of the vector (1, 2, 3).                   */
int main(void)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ThreeVector p;
    double x, y, z, norm;

    /*  Set the x, y, and z values to 1, 2, and 3, respectively.              */
    x = 1.0;
    y = 2.0;
    z = 3.0;

    /*  Set p to the vector (x, y, z).                                        */
    p = rssringoccs_ThreeVector_Rect(x, y, z);

    /*  Compute the norm of p.                                                */
    norm = rssringoccs_ThreeVector_Euclidean_Norm(p);

    /*  Print the result:                                                     */
    printf("||(%f, %f, %f)|| = %f\n", x, y, z, norm);

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc euclidean_norm_3d_example.c -o test -lrssringoccs                 *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), change *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              euclidean_norm_3d_example.c -o test -lrssringoccs             *
 *                                                                            *
 *  This example is also C89 compliant and compiles with the following flags: *
 *                                                                            *
 *      gcc -Wconversion -pedantic -Wall -Wextra -std=c89 -ansi               *
 *          -Wpedantic euclidean_norm_3d_example.c -o test -lrssringoccs      *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      ||(1.000000, 2.000000, 3.000000)|| = 3.741657                         *
 *  Which is equal to the square root of 14.                                  *
 ******************************************************************************/
