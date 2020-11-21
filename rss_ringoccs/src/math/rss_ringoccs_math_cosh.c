/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Cosh_Float(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float cosh_x;

    /*  The definition of cosh(x) is [exp(x) + exp(-x)]/2, so return this.    */
    cosh_x = 0.5*(rssringoccs_Exp_Float(x) + rssringoccs_Exp_Float(-x));
    return cosh_x;
}

double rssringoccs_Cosh_Double(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double cosh_x;

    /*  The definition of sinh(x) is [exp(x) + exp(-x)]/2, so return this.    */
    cosh_x = 0.5*(rssringoccs_Exp_Double(x) + rssringoccs_Exp_Double(-x));
    return cosh_x;
}

long double rssringoccs_Cosh_LongDouble(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double cosh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    cosh_x = 0.5*(rssringoccs_Exp_LongDouble(x)+rssringoccs_Exp_LongDouble(-x));
    return cosh_x;
}
