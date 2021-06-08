#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double
rssringoccs_Fresnel_Scale(double lambda, double d, double phi, double b)
{
    double cb_2_sp_2, sb_2, f_scale;

    cb_2_sp_2  = cos(b)*sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2 = cos(b);
    sb_2 *= sb_2;

    f_scale = sqrt(0.5*lambda*d*(1.0 - cb_2_sp_2)/sb_2);
    return f_scale;
}
