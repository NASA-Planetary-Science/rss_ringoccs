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
 *      Provide tests for the accuracy and efficiency of rss_ringoccs         *
 *      sine function compared to the one provided in math.h.                 *
 *  NOTE:                                                                     *
 *      If rss_ringoccs was built with the macro                              *
 *      __RSS_RINGOCCS_USE_TRIG_ALGORITHMS__ was set to 0 in                  *
 *      rss_ringoccs_config.h, then this test is redundant since sinl and     *
 *      rssringoccs_LDouble_Sin are the same thing.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 12, 2020                                             *
 ******************************************************************************/

/*  The sine functions are found here.                                        */
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <math.h>

/*  The comparison functions are found here.                                  *
 *  NOTE:                                                                     *
 *      You will need to compile rss_ringoccs_compare_funcs.c. This file and  *
 *      the functions found in rss_ringoccs_compare_funcs.h are NOT found in  *
 *      librssringoccs. We can compile via:                                   *
 *                                                                            *
 *          gcc -O3 -pedantic -Wall -Wconversion -Wextra -Wpedantic           *
 *              rss_ringoccs_compare_funcs.c -shared                          *
 *                  -o librssringoccs_compare.so                              *
 *                                                                            *
 *      In the examples below we placed the output file in /usr/local/lib/:   *
 *                                                                            *
 *          mv librssringoccs_compare.so /usr/local/lib/                      *
 *                                                                            *
 *      We can then link via -lrssringoccs_compare (see below).               */
#include "../rss_ringoccs_compare_funcs.h"

/*  Routine for comparing sinl with rssringoccs_LDouble_Sin.                  */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    long double start = -100.0L;
    long double end   =  100.0L;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e8;

    /*  Use the compare function to test rssringoccs_LDouble_Sin against sinl.*/
    rssringoccs_Compare_LDouble_Funcs("C99", sinl,
                                      "rss_ringoccs", rssringoccs_LDouble_Sin,
                                      start, end, N);

    return 0;
}
/*  End of main.                                                              */

