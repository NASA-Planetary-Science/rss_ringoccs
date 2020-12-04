/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Float_Cosh(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float cosh_x;

    /*  The definition of cosh(x) is [exp(x) + exp(-x)]/2, so return this.    */
    cosh_x = 0.5F*(rssringoccs_Float_Exp(x) + rssringoccs_Float_Exp(-x));
    return cosh_x;
}

double rssringoccs_Double_Cosh(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double cosh_x;

    /*  The definition of sinh(x) is [exp(x) + exp(-x)]/2, so return this.    */
    cosh_x = 0.5*(rssringoccs_Double_Exp(x) + rssringoccs_Double_Exp(-x));
    return cosh_x;
}

long double rssringoccs_LDouble_Cosh(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double cosh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    cosh_x = 0.5L*(rssringoccs_LDouble_Exp(x)+rssringoccs_LDouble_Exp(-x));
    return cosh_x;
}
