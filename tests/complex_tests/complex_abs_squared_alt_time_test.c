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
 *  Date:       December 20, 2020                                             *
 ******************************************************************************/

/*  All necessary complex data types found in these two files.                */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <complex.h>

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

/*  C99 does not provide an abs squared function, so let's create one using   *
 *  the built-in _Complex data type.                                          */
static double cabs_sq(_Complex double z)
{
    /*  Declare necessary variables.                                          */
    double x, y, abs_sq;

    /*  Use the creal and cimag functions found in complex.h to extract the   *
     *  real and imaginary parts from the input z.                            */
    x = creal(z);
    y = cimag(z);

    /*  |z|^2 = x^2 + y^2 so return this.                                     */
    abs_sq = x*x + y*y;
    return abs_sq;
}

/*  The rssringoccs_CDouble_Abs_Squared function was written with clarity in  *
 *  mind. We can save calls to the functions rssringoccs_CDouble_Real_Part    *
 *  and rssringoccs_CDouble_Imag_Part by directly accessing the data member   *
 *  of the rssringoccs_ComplexDouble struct. This looks a little more cryptic *
 *  and is most likely an example of premature optimization. Nonetheless, the *
 *  time test is provided here.                                               */
static double rss_cabs_sq(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables.                                          */
    double x, y, abs_sq;

    /*  Access the real and imaginary parts of the struct z directly, rather  *
     *  than making calls to rssringoccs_CDouble_Real/Imag_Part.              */
    x = z.dat[0];
    y = z.dat[1];

    /*  |z|^2 = x^2 + y^2 so return this.                                     */
    abs_sq = x*x + y*y;
    return abs_sq;
}

/*  Routine for comparing cabs_sq with rss_cabs_sq.                           */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    double start = -1.0;
    double end   =  1.0;

    /*  We'll test on a square grid of 100 million points from (start, start) *
     *  the (end, end) in the complex plane.                                  */
    unsigned long N = 1e4;

    /*  Use the compare function found in rss_ringoccs_compare_funcs.h.       */
    rssringoccs_Compare_Real_CDouble_Funcs("rss_ringoccs", rss_cabs_sq,
                                           "C99", cabs_sq, start, end, N);

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  Compileable with:                                                         *
 *      gcc -O3 -Wall -Wpedantic -Wextra -pedantic -pedantic-errors           *
 *          -std=c99 complex_abs_squared_alt_time_test.c -o test              *
 *              -lrssringoccs -lrssringoccs_compare                           *
 *  Output (iMac 2017 3.4 GHz Intel Quad-Core i5):                            *
 *      rss_ringoccs: 0.692513                                                *
 *      C99: 0.664646                                                         *
 *      Max Error: 0.0000000000000000                                         *
 *  With -O3 optimization:                                                    *
 *      rss_ringoccs: 0.633626                                                *
 *      C99: 0.560142                                                         *
 *      Max Error: 0.0000000000000000                                         *
 ******************************************************************************/

/*  For comparison, the times for rssringoccs_CDouble_Abs_Squared are:        *
 *      Optimization:       0.816004                                          *
 *      w/o optimization:   0.817405                                          *
 *  These tests were performed on a grid of 100 million points in the complex *
 *  plane and took less than a second. Avoiding calls to the two functions    *
 *  rssringoccs_CDouble_Real/Imag_Part saves 0.2 seconds over 100 million     *
 *  computations, which is barely noticeable.                                 */
