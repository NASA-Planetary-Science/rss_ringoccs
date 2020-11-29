#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float Fresnel_Scale_Float(float lambda, float d, float phi, float b)
{
    float cb_2_sp_2, sb_2;

    cb_2_sp_2  = cosf(b)*sinf(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sinf(b);
    sb_2 *= sb_2;

    return sqrtf(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}

extern double
Fresnel_Scale_Double(double lambda, double d, double phi, double b)
{
    double cb_2_sp_2, sb_2;

    cb_2_sp_2  = cos(b)*sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sin(b);
    sb_2 *= sb_2;

    return sqrt(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}


extern long double
Fresnel_Scale_LongDouble(long double lambda, long double d,
                          long double phi, long double b)
{
    long double cb_2_sp_2, sb_2;

    cb_2_sp_2  = cosl(b)*sinl(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sinl(b);
    sb_2 *= sb_2;

    return sqrtl(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}
