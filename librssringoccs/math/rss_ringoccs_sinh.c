/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Float_Sinh(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5F*(rssringoccs_Float_Exp(x) - rssringoccs_Float_Exp(-x));
    return sinh_x;
}

double rssringoccs_Double_Sinh(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5*(rssringoccs_Double_Exp(x) - rssringoccs_Double_Exp(-x));
    return sinh_x;
}

long double rssringoccs_LDouble_Sinh(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5L*(rssringoccs_LDouble_Exp(x)-rssringoccs_LDouble_Exp(-x));
    return sinh_x;
}
