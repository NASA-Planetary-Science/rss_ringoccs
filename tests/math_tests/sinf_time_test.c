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
 *      rss_ringoccs_config.h, then this test is redundant since sinf and     *
 *      rssringoccs_Float_Sin are the same thing.                             *
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

/*  Routine for comparing sinf with rssringoccs_Float_Sin.                    */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    float start = -100.0F;
    float end   =  100.0F;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e8;

    /*  Use the compare function to test rssringoccs_Float_Sin against sinf.  */
    rssringoccs_Compare_Float_Funcs("C99", sinf,
                                    "rss_ringoccs", rssringoccs_Float_Sin,
                                    start, end, N);

    return 0;
}
/*  End of main.                                                              */

/*  This was compiled with various options on an iMac 2017 3.4GHz quad-core   *
 *  running MacOS Catalina 10.15.7. It produced the following times:          *
 *      NOTE: On MacOS gcc is aliased to LLVM's clang:                        *
 *          gcc --version                                                     *
 *          Apple clang version 12.0.0 (clang-1200.0.32.27)                   *
 *      This is NOT the regular gcc from GNU. To use gcc on apple devices     *
 *      requires homebrew.                                                    *
 *  c99 option, -O3 opimization:                                              *
 *      gcc -pedantic -Wall -Wextra -Wpedantic -Wconversion -std=c99 -O3      *
 *              sinf_time_test.c -o test -lrssringoccs -lrssringoccs_compare  *
 *      C99:          0.899424                                                *
 *      rss_ringoccs: 1.315589                                                *
 *      Max Error: 0.00000000                                                 *
 *  c99 option, no optimization.                                              *
 *      gcc -pedantic -Wall -Wextra -Wpedantic -Wconversion -std=c99          *
 *              sinf_time_test.c -o test -lrssringoccs -lrssringoccs_compare  *
 *      C99:          0.877683                                                *
 *      rss_ringoccs: 1.336257                                                *
 *      Max Error: 0.00000000                                                 *
 *  NOTE:                                                                     *
 *      These times will differ on different devices and on different         *
 *      compilers.                                                            */
