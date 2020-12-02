/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  The "double" version of cos is defined in both C89 and C99 math.h so we   *
 *  only need to alias this function.                                         */
double rssringoccs_Double_Abs(double x)
{
    return fabs(x);
}

#if __HAS_C99_MATH_H__ == 0

/*  C89 math.h does not have cosf or cosfl, so we'll need to provide these to  *
 *  make the code forward compatible. We'll do this in a very simple manner.  */
float rssringoccs_Float_Abs(float x)
{
    return fabs((double)x);
}

long double rssringoccs_LDouble_Abs(long double x)
{
    return fabs((double)x);
}

/*  Now have the functions declared in rss_ringoccs_math.h point to these.    */
#else

/*  C99 provides float and long double support for their math functions, so   *
 *  simply use to these.                                                      */
float rssringoccs_Float_Abs(float x)
{
    return fabsf(x);
}

long double rssringoccs_LDouble_Abs(long double x)
{
    return fabsl(x);
}
#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */
