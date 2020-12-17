

/*  The absolute value functions are found here.                              */
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
    double start = -1.0;
    double end   =  1.0;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e4;

    /*  Use the compare function to test rssringoccs_Double_Abs against fabs. */
    rssringoccs_Compare_CFloat_Funcs("rss_ringoccs", rssringoccs_CFloat_Cos,
                                     "C99", ccosf, start, end, N);

    return 0;
}