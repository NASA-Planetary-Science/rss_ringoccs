/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  The "double" version of cos is defined in both C89 and C99 math.h so we   *
 *  only need to alias this function.                                         */
double rssringoccs_Double_Arctan2(double y, double x)
{
    return atan2(y, x);
}

#if __HAS_C99_MATH_H__ == 0

/*  C89 math.h does not have cosf or cosfl, so we'll need to provide these to  *
 *  make the code forward compatible. We'll do this in a very simple manner.  */
float rssringoccs_Float_Arctan2(float y, float x)
{
    return atan2((double)y, (double)x);
}

long double rssringoccs_LDouble_Arctan2(long double y, long double x)
{
    return atan2((double)y, (double)x);
}

/*  Now have the functions declared in rss_ringoccs_math.h point to these.    */
#else

/*  C99 provides float and long double support for their math functions, so   *
 *  simply use to these.                                                      */
float rssringoccs_Float_Arctan2(float y, float x)
{
    return atan2f(y, x);
}

long double rssringoccs_LDouble_Arctan2(long double y, long double x)
{
    return atan2l(y, x);
}

#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */
