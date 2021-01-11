/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Float_Tanh(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_Float_Sinh(x) / rssringoccs_Float_Cosh(-x);
    return tanh_x;
}

double rssringoccs_Double_Tanh(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_Double_Sinh(x) / rssringoccs_Double_Cosh(-x);
    return tanh_x;
}

long double rssringoccs_LDouble_Tanh(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_LDouble_Sinh(x) / rssringoccs_LDouble_Cosh(-x);
    return tanh_x;
}
