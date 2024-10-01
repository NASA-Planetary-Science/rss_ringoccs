#include <libtmpl/include/tmpl_vec3.h>


#include <rss_ringoccs/include/rss_ringoccs_cspice.h>

#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>

#define RSSRINGOCCS_RADIANS_PER_DEGREE (+5.72957795130823208767981548141051E+01)
#define RSSRINGOCCS_PI_BY_TWO (+1.570796326794896619231321691639751442098E+00)

double
rssringoccs_Double_Ring_Opening_Angle(const double time,
                                      const tmpl_ThreeVectorDouble * const nhat,
                                      const char * const spacecraft,
                                      const char * const dsn,
                                      const char * const kernels,
                                      const char * const ref)
{
    double ltime, opening_radians;
    tmpl_ThreeVectorDouble ptarg;

    if (kernels)
    {
        kclear_c();
        furnsh_c(kernels);
    }

    spkpos_c(dsn, time, ref, "CN", spacecraft, ptarg.dat, &ltime);

    opening_radians = RSSRINGOCCS_PI_BY_TWO - tmpl_3DDouble_Angle(&ptarg, nhat);
    return opening_radians * RSSRINGOCCS_RADIANS_PER_DEGREE;
}

/*  Undefine everything in case someone wants to #include this file.          */
#undef RSSRINGOCCS_PI_BY_TWO
#undef RSSRINGOCCS_RADIANS_PER_DEGREE
