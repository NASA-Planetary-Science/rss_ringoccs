/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Float_Erfc(float x)
{
    float erfc;
    erfc = 1.0 - rssringoccs_Float_Erf(x);

    return erfc;
}

double rssringoccs_Double_Erfc(double x)
{
    double erfc;
    erfc = 1.0 - rssringoccs_Double_Erf(x);

    return erfc;
}

long double rssringoccs_LongDouble_Erfc(long double x)
{
    long double erfc;
    erfc = 1.0 - rssringoccs_LongDouble_Erf(x);

    return erfc;
}
