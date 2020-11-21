/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Tanh_Float(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_Sinh_Float(x) / rssringoccs_Cosh_Float(-x);
    return tanh_x;
}

double rssringoccs_Tanh_Double(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_Sinh_Double(x) / rssringoccs_Cosh_Double(-x);
    return tanh_x;
}

long double rssringoccs_Tanh_LongDouble(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float tanh_x;

    /*  The definition of than(x) is sinh(x)/cosh(x), so return this.         */
    tanh_x = rssringoccs_Sinh_LongDouble(x) / rssringoccs_Cosh_LongDouble(-x);
    return tanh_x;
}
