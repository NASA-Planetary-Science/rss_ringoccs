/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

double rssringoccs_Double_Tan(double x)
{
    return tan(x);
}

#if __HAS_C99_MATH_H__ == 0

float rssringoccs_Float_Tan(float x)
{
    return (float)tan((double)x);
}

long double rssringoccs_LDouble_Tan(long double x)
{
    return (long double)tan((double)x);
}

/*  Now have the functions declared in rss_ringoccs_math.h point to these.    */
#else

float rssringoccs_Float_Tan(float x)
{
    return tanf(x);
}

long double rssringoccs_LDouble_Tan(long double x)
{
    return tanl(x);
}
#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */
