#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>
#include <libtmpl/include/tmpl.h>
#include <cspice/include/SpiceUsr.h>

double
rssringoccs_Double_Ring_Opening_Angle(const double time,
                                      const tmpl_ThreeVectorDouble * const nhat,
                                      const char *const spacecraft,
                                      const char *const dsn,
                                      const char *const kernels,
                                      const char *const ref)
{
    double ltime, opening_radians;
    tmpl_ThreeVectorDouble ptarg;

    if (kernels)
    {
        kclear_c();
        furnsh_c(kernels);
    }

    spkpos_c(dsn, time, ref, "CN", spacecraft, ptarg.dat, &ltime);

    opening_radians = tmpl_Pi_By_Two - vsep_c(ptarg.dat, nhat->dat);
    return opening_radians * tmpl_Rad_to_Deg;
}
