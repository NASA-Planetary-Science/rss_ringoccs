#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float
rssringoccs_Float_Fresnel_Scale(float lambda, float d, float phi, float b)
{
    float cb_2_sp_2, sb_2, f_scale;

    cb_2_sp_2  = rssringoccs_Float_Cos(b)*rssringoccs_Float_Sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = rssringoccs_Float_Sin(b);
    sb_2 *= sb_2;

    f_scale = rssringoccs_Float_Sqrt(0.5F*lambda*d*(1.0F - cb_2_sp_2)/sb_2);
    return f_scale;
}

extern double
rssringoccs_Double_Fresnel_Scale(double lambda, double d, double phi, double b)
{
    double cb_2_sp_2, sb_2, f_scale;

    cb_2_sp_2  = rssringoccs_Double_Cos(b)*rssringoccs_Double_Sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = rssringoccs_Double_Sin(b);
    sb_2 *= sb_2;

    f_scale = rssringoccs_Double_Sqrt(0.5*lambda*d*(1.0 - cb_2_sp_2)/sb_2);
    return f_scale;
}


extern long double
rssringoccs_LDouble_Fresnel_Scale(long double lambda, long double d,
                                  long double phi, long double b)
{
    long double cb_2_sp_2, sb_2, f_scale;

    cb_2_sp_2  = rssringoccs_LDouble_Cos(b)*rssringoccs_LDouble_Sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = rssringoccs_LDouble_Sin(b);
    sb_2 *= sb_2;

    f_scale = rssringoccs_LDouble_Sqrt(0.5L*lambda*d*(1.0L - cb_2_sp_2)/sb_2);
    return f_scale;
}
