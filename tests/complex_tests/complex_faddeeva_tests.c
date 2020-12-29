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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 19, 2020                                             *
 ******************************************************************************/

/******************************************************************************
 *  NOTES:                                                                    *
 *      The error functions grows faster than exp(y^2) for z = 0 + iy.        *
 *      Because of this, the test here calculates the relative error, which   *
 *      0.0000000000000004 (see below). The worst absolute error is           *
 *      618970019642690137449562112.0000000000000000. Since the value of this *
 *      function has magnitude ~10^42 at this point, we see why the absolute  *
 *      error is so huge.                                                     *
 ******************************************************************************/

/*  rss_ringoccs complex routines found here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Complex error function provided by the MIT Faddeeva package found here.   *
 *  You must have the library installed and in your path to use this test.    */
#include <cerf.h>

/*  The comparison functions are found here.                                  *
 *  NOTE:                                                                     *
 *      You will need to compile rss_ringoccs_compare_funcs.c. This file and  *
 *      the functions found in rss_ringoccs_compare_funcs.h are NOT found in  *
 *      librssringoccs. We can compile via:                                   *
 *                                                                            *
 *          gcc -O3 -pedantic -Wall -Wconversion -Wextra -Wpedantic           *
 *              rss_ringoccs_compare_funcs.c -shared                          *
 *                  -o librssringoccs_compare.so -lrssringoccs                *
 *                                                                            *
 *      You will need to have librssringoccs already built so we can link it  *
 *      to this new comparison library. In the examples below we placed the   *
 *      output file in /usr/local/lib/:                                       *
 *                                                                            *
 *          mv librssringoccs_compare.so /usr/local/lib/                      *
 *                                                                            *
 *      You may need to add sudo to move files in /usr/.                      *
 *      We can then link via -lrssringoccs_compare (see below).               */
#include "../rss_ringoccs_compare_funcs.h"

/*  Routine for testing rssringoccs_CDouble_Erf.                              */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    double start = -10.0;
    double end   =  10.0;

    /*  We'll test on a square grid of 100 million points from (start, start) *
     *  the (end, end) in the complex plane.                                  */
    unsigned long N = 1e3;

    /*  Use the compare function found in rss_ringoccs_compare_funcs.h.       */
    rssringoccs_RelCompare_CDouble_Funcs("rss_ringoccs",
                                         rssringoccs_CDouble_Faddeeva,
                                         "libcerf", w_of_z, start, end, N);



    /*  And run an accuracy test over a larger region of the complex plane.   */
    start = -5000.0;
    end   =  5000.0;
    N = 1e3;

    rssringoccs_Accuracy_CDouble_Funcs("rss_ringoccs",
                                        rssringoccs_CDouble_Faddeeva,
                                        "libcerf", w_of_z, start, end, N);
    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  Compileable with:                                                         *
 *      gcc -O3 -Wall -Wpedantic -Wextra -pedantic -pedantic-errors           *
 *          -std=c99 complex_erf_time_test.c -o test -lrssringoccs            *
 *              -lrssringoccs_compare -lcerf                                  *
 *      Don't forget to link libcerf as well!                                 *
 *  Output (iMac 2017 3.4 GHz Intel Quad-Core i5):                            *
 *      rss_ringoccs: 1.922353                                                *
 *      libcerf: 1.076348                                                     *
 *      Max Relative Error: 0.0000000000000004                                *
 *  With -O3 optimization:                                                    *
 *      rss_ringoccs: 1.928797                                                *
 *      libcerf: 1.182088                                                     *
 *      Max Relative Error: 0.0000000000000004                                *
 ******************************************************************************/
