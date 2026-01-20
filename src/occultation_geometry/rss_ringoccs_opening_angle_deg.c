#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <rss_ringoccs/include/rss_ringoccs_cspice.h>
#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>

double
rssringoccs_Double_Ring_Opening_Angle_Deg(
    const double time,
    const tmpl_ThreeVectorDouble * TMPL_RESTRICT const nhat,
    const char * TMPL_RESTRICT const spacecraft,
    const char * TMPL_RESTRICT const dsn,
    const char * TMPL_RESTRICT const ref
)
{
    const char * const abcorr = "CN";
    double ltime, opening_angle;
    tmpl_ThreeVectorDouble ptarg;

    spkpos_c(dsn, time, ref, abcorr, spacecraft, ptarg.dat, &ltime);

    opening_angle = tmpl_double_pi_by_two - tmpl_3DDouble_Angle(&ptarg, nhat);
    return opening_angle * tmpl_double_rad_to_deg;
}
