/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include "rss_ringoccs_math.h"

float rssringoccs_Sinh_Float(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5*(rssringoccs_Exp_Float(x) - rssringoccs_Exp_Float(-x));
    return sinh_x;
}

double rssringoccs_Sinh_Double(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5*(rssringoccs_Exp_Double(x) - rssringoccs_Exp_Double(-x));
    return sinh_x;
}

long double rssringoccs_Sinh_LongDouble(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double sinh_x;

    /*  The definition of sinh(x) is [exp(x) - exp(-x)]/2, so return this.    */
    sinh_x = 0.5*(rssringoccs_Exp_LongDouble(x)-rssringoccs_Exp_LongDouble(-x));
    return sinh_x;
}

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
