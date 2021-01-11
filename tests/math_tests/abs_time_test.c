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
 *      absolute value function compared to the one provided by the C99       *
 *      standard in math.h.                                                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  The absolute value functions are found here.                              */
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

/*  Routine for comparing fabs with rssringoccs_Double_Abs.                   */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    double start = -100.0;
    double end   =  100.0;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e8;

    /*  Use the compare function to test rssringoccs_Double_Abs against fabs. */
    rssringoccs_Compare_Double_Funcs("C99", fabs,
                                     "rss_ringoccs", rssringoccs_Double_Abs,
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
 *  c89 option, -O3 opimization:                                              *
 *      gcc -Wall -Wextra -Wpedantic -Wconversion -std=c89 -ansi -O3          *
 *              abs_time_test.c -o test -lrssringoccs -lrssringoccs_compare   *
 *      C99:          0.562831                                                *
 *      rss_ringoccs: 0.581634                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c89 option, no optimization.                                              *
 *      gcc -pedantic -Wall -Wextra -Wpedantic -Wconversion -std=c89 -ansi    *
 *              abs_time_test.c -o test -lrssringoccs -lrssringoccs_compare   *
 *      C99:          0.563751                                                *
 *      rss_ringoccs: 0.612806                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  NOTE:                                                                     *
 *      These times will differ on different devices and on different         *
 *      compilers.                                                            */
