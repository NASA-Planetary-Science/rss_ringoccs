/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/src/math/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include "rss_ringoccs_complex.h"

static unsigned int __cerf_taylor_deg = 34;
static double __cerf_taylor_coeffs[35] = {
    1.0,                    -0.3333333333333333,
    0.1,                    -0.023809523809523808,
    0.004629629629629629,   -0.0007575757575757576,
    0.00010683760683760684, -1.3227513227513228e-05,
    1.4589169000933706e-06, -1.4503852223150468e-07,
    1.3122532963802806e-08, -1.0892221037148573e-09,
    8.35070279514724e-11,   -5.9477940136376354e-12,
    3.9554295164585257e-13, -2.466827010264457e-14,
    1.4483264643598138e-15, -8.032735012415773e-17,
    4.221407288807088e-18,  -2.107855191442136e-19,
    1.0025164934907719e-20, -4.5518467589282e-22,
    1.977064753877905e-23,  -8.230149299214221e-25,
    3.289260349175752e-26,  -1.2641078988989164e-27,
    4.6784835155184856e-29, -1.669761793417372e-30,
    5.754191643982172e-32,  -1.9169428621097826e-33,
    6.180307588222796e-35,  -1.930357208815108e-36,
    5.846755007468836e-38,  -1.7188560628017835e-39,
    4.908923964523423e-41
};

static unsigned int __cerf_asym_deg = 5;
static double __cerf_asym_coeffs[6] = {
    1.0,    -0.5,
    0.75,   -1.875,
    6.5625, -29.5313
};

rssringoccs_ComplexDouble
rssringoccs_Complex_Erf(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double abs_z, real;
    rssringoccs_ComplexDouble erfz, arg, scale;

    abs_z = rssringoccs_Complex_Abs(z);
    arg = rssringoccs_Complex_Multiply(z, z);

    if (abs_z < 3.0)
    {
        erfz = rssringoccs_Complex_Poly_Real_Coeffs(
            __cerf_taylor_coeffs, __cerf_taylor_deg, arg
        );

        erfz = rssringoccs_Complex_Multiply(z, erfz);
        erfz = rssringoccs_Complex_Scale(TWO_BY_SQRT_PI, erfz);
    }
    else
    {
        real = rssringoccs_Complex_Real_Part(z);
        if (real < 0.0)
            z = rssringoccs_Complex_Scale(-1.0, z);

        scale = rssringoccs_Complex_Scale(-1.0, arg);
        scale = rssringoccs_Complex_Exp(scale);
        scale = rssringoccs_Complex_Scale(ONE_BY_SQRT_PI, scale);
        arg = rssringoccs_Complex_Reciprocal(arg);

        erfz = rssringoccs_Complex_Poly_Real_Coeffs(
            __cerf_asym_coeffs, __cerf_asym_deg, arg
        );

        erfz = rssringoccs_Complex_Divide(erfz, z);
        erfz = rssringoccs_Complex_Multiply(erfz, scale);
        erfz = rssringoccs_Complex_Subtract(rssringoccs_Complex_One, erfz);

        if (real < 0.0)
            erfz = rssringoccs_Complex_Scale(-1.0, erfz);
    }
    return erfz;
}
