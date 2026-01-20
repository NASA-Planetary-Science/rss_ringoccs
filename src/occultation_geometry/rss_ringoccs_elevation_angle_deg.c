#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <rss_ringoccs/include/rss_ringoccs_cspice.h>
#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>

/*  TODO:
 *      Investigate need to swap Malargue's number with name:
 *
 *      if observer == '398958':
 *          observer = 'MALARGUE'
 */

double
rssringoccs_Double_Elevation_Angle_Deg(
    const double time,
    const char * TMPL_RESTRICT const target,
    const char * TMPL_RESTRICT const observer,
    const char * TMPL_RESTRICT const ref
)
{
    const char * const abcorr = "CN";
    const char * const planet = "EARTH";
    double ltime1, ltime2;
    double elevation_angle, separation_angle;
    tmpl_ThreeVectorDouble ptarg1, ptarg2;

    spkpos_c(target, time, ref, abcorr, observer, ptarg1.dat, &ltime1);
    spkpos_c(observer, time, ref, abcorr, planet, ptarg2.dat, &ltime2);

    separation_angle = tmpl_3DDouble_Angle(&ptarg1, &ptarg2);
    elevation_angle = tmpl_double_pi_by_two - separation_angle;
    return elevation_angle * tmpl_double_rad_to_deg;
}
