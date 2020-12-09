/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include "rss_ringoccs_math_private.h"

#if __RSS_RINGOCCS_USE_TRIG_ALGORITHMS__ != 0

float rssringoccs_Float_Sin(float x)
{
    float arg, abs_x, sgn_x, cx, cdx, sx, sdx, dx;
    unsigned int arg_100_int;

    if (x >= 0.0F)
    {
        abs_x = x;
        sgn_x = 1.0F;
    }
    else
    {
        abs_x = -x;
        sgn_x = -1.0F;
    }

    arg = (float)fmod((double)abs_x, rssringoccs_Two_Pi);

    if (arg > rssringoccs_One_Pi_F)
    {
        sgn_x *= -1.0F;
        arg -= rssringoccs_One_Pi_F;
    }

    if (arg < rssringoccs_Pi_By_Four_F)
    {
        arg_100_int = (unsigned int)(100.0F*arg);
        dx = arg - 0.01F*arg_100_int;

        sx  = rssringoccs_sinf_table(arg_100_int);
        cdx = rssringoccs_do_cosf(dx);

        cx  = rssringoccs_cosf_table(arg_100_int);
        sdx = rssringoccs_do_sinf(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
    else if (arg < rssringoccs_Pi_By_Two_F)
    {
        arg = rssringoccs_Pi_By_Two_F - arg;
        arg_100_int = (unsigned int)(100.0F*arg);
        dx = arg - 0.01F*arg_100_int;

        cx  = rssringoccs_cosf_table(arg_100_int);
        cdx = rssringoccs_do_cosf(dx);

        sx  = rssringoccs_sinf_table(arg_100_int);
        sdx = rssringoccs_do_sinf(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else if (arg < rssringoccs_Three_Pi_By_Four_F)
    {
        arg -= rssringoccs_Pi_By_Two_F;
        arg_100_int = (unsigned int)(100.0F*arg);
        dx = arg - 0.01F*arg_100_int;

        cx  = rssringoccs_cosf_table(arg_100_int);
        cdx = rssringoccs_do_cosf(dx);

        sx  = rssringoccs_sinf_table(arg_100_int);
        sdx = rssringoccs_do_sinf(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else
    {
        arg = rssringoccs_One_Pi_F - arg;
        arg_100_int = (unsigned int)(100.0F*arg);
        dx = arg - 0.01F*arg_100_int;

        sx  = rssringoccs_sinf_table(arg_100_int);
        cdx = rssringoccs_do_cosf(dx);

        cx  = rssringoccs_cosf_table(arg_100_int);
        sdx = rssringoccs_do_sinf(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
}

double rssringoccs_Double_Sin(double x)
{
    double arg, abs_x, sgn_x, cx, cdx, sx, sdx, dx;
    unsigned int arg_100_int;

    if (x >= 0.0)
    {
        abs_x = x;
        sgn_x = 1.0;
    }
    else
    {
        abs_x = -x;
        sgn_x = -1.0;
    }

    arg = fmod(abs_x, rssringoccs_Two_Pi);

    if (arg > rssringoccs_One_Pi)
    {
        sgn_x *= -1.0;
        arg -= rssringoccs_One_Pi;
    }

    if (arg < rssringoccs_Pi_By_Four)
    {
        arg_100_int = (unsigned int)(100.0*arg);
        dx = arg - 0.01*arg_100_int;

        sx  = rssringoccs_sin_table(arg_100_int);
        cdx = rssringoccs_do_cos(dx);

        cx  = rssringoccs_cos_table(arg_100_int);
        sdx = rssringoccs_do_sin(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
    else if (arg < rssringoccs_Pi_By_Two)
    {
        arg = rssringoccs_Pi_By_Two - arg;
        arg_100_int = (unsigned int)(100.0*arg);
        dx = arg - 0.01*arg_100_int;

        cx  = rssringoccs_cos_table(arg_100_int);
        cdx = rssringoccs_do_cos(dx);

        sx  = rssringoccs_sin_table(arg_100_int);
        sdx = rssringoccs_do_sin(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else if (arg < rssringoccs_Three_Pi_By_Four)
    {
        arg -= rssringoccs_Pi_By_Two;
        arg_100_int = (unsigned int)(100.0*arg);
        dx = arg - 0.01*arg_100_int;

        cx  = rssringoccs_cos_table(arg_100_int);
        cdx = rssringoccs_do_cos(dx);

        sx  = rssringoccs_sin_table(arg_100_int);
        sdx = rssringoccs_do_sin(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else
    {
        arg = rssringoccs_One_Pi - arg;
        arg_100_int = (unsigned int)(100.0*arg);
        dx = arg - 0.01*arg_100_int;

        sx  = rssringoccs_sin_table(arg_100_int);
        cdx = rssringoccs_do_cos(dx);

        cx  = rssringoccs_cos_table(arg_100_int);
        sdx = rssringoccs_do_sin(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
}

long double rssringoccs_LDouble_Sin(long double x)
{
    long double arg, abs_x, sgn_x, cx, cdx, sx, sdx, dx;
    unsigned int arg_100_int;

    if (x >= 0.0L)
    {
        abs_x = x;
        sgn_x = 1.0L;
    }
    else
    {
        abs_x = -x;
        sgn_x = -1.0L;
    }

    if (abs_x > rssringoccs_Two_Pi_L)
        arg = (long double)fmod((double)abs_x, rssringoccs_Two_Pi);
    else
        arg = abs_x;

    if (arg > rssringoccs_One_Pi_L)
    {
        sgn_x *= -1.0L;
        arg -= rssringoccs_One_Pi_L;
    }

    if (arg < rssringoccs_Pi_By_Four_L)
    {
        arg_100_int = (unsigned int)(100.0L*arg);
        dx = arg - 0.01L*arg_100_int;

        sx  = rssringoccs_sinl_table(arg_100_int);
        cdx = rssringoccs_do_cosl(dx);

        cx  = rssringoccs_cosl_table(arg_100_int);
        sdx = rssringoccs_do_sinl(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
    else if (arg < rssringoccs_Pi_By_Two_L)
    {
        arg = rssringoccs_Pi_By_Two_L - arg;
        arg_100_int = (unsigned int)(100.0L*arg);
        dx = arg - 0.01L*arg_100_int;

        cx  = rssringoccs_cosl_table(arg_100_int);
        cdx = rssringoccs_do_cosl(dx);

        sx  = rssringoccs_sinl_table(arg_100_int);
        sdx = rssringoccs_do_sinl(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else if (arg < rssringoccs_Three_Pi_By_Four_L)
    {
        arg -= rssringoccs_Pi_By_Two_L;
        arg_100_int = (unsigned int)(100.0L*arg);
        dx = arg - 0.01L*arg_100_int;

        cx  = rssringoccs_cosl_table(arg_100_int);
        cdx = rssringoccs_do_cosl(dx);

        sx  = rssringoccs_sinl_table(arg_100_int);
        sdx = rssringoccs_do_sinl(dx);

        return sgn_x*(cx*cdx - sx*sdx);
    }
    else
    {
        arg = rssringoccs_One_Pi_L - arg;
        arg_100_int = (unsigned int)(100.0L*arg);
        dx = arg - 0.01L*arg_100_int;

        sx  = rssringoccs_sinl_table(arg_100_int);
        cdx = rssringoccs_do_cosl(dx);

        cx  = rssringoccs_cosl_table(arg_100_int);
        sdx = rssringoccs_do_sinl(dx);

        return sgn_x*(sx*cdx + cx*sdx);
    }
}

#else

/*  The "double" version of sin is defined in both C89 and C99 math.h so we   *
 *  only need to alias this function.                                         */

/*  Double precision sine function (sin equivalent).                          */
double rssringoccs_Double_Sin(double x)
{
    return sin(x);
}
/*  End of rssringoccs_Double_Sin.                                            */


#if __HAS_C99_MATH_H__ == 0

/*  Single precision sine function (sinf equivalent).                         */
float rssringoccs_Float_Sin(float x)
{
    return (float)sin((double)x);
}
/*  End of rssringoccs_Float_Sin.                                             */

/*  Long double precision sine function (sinl equivalent).                    */
long double rssringoccs_LDouble_Sin(long double x)
{
    return (long double)sin((double)x);
}
/*  End of rssringoccs_LDouble_Sin.                                           */

#else

/*  C99 provides float and long double support for their math functions, so   *
 *  simply use to these.                                                      */
float rssringoccs_Float_Sin(float x)
{
    return sinf(x);
}

long double rssringoccs_LDouble_Sin(long double x)
{
    return sinl(x);
}
#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */

#endif
/*  End of #if __RSS_RINGOCCS_USE_TRIG_ALGORITHMS__ != 0.                     */

