#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_compare_funcs.h>

int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    float start = 0.0F;
    float end   = 100.0F;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e8;

    /*  Use the compare function to test rssringoccs_Double_Abs against fabs. */
    rssringoccs_Compare_Float_Funcs("C99", logf,
                                    "rss_ringoccs", rssringoccs_Float_Log,
                                    start, end, N);

    return 0;
}
