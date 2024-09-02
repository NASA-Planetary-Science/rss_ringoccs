#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>
#include <libtmpl/include/tmpl.h>

extern void kclear_c(void);
extern void furnsh_c(const char *kernels);
extern void
spkpos_c(const char *targ, double et, const char *ref, const char *abcorr,
         const char *obs, double ptarg[3], double *lt);

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

    opening_radians = tmpl_Pi_By_Two - tmpl_3DDouble_Angle(&ptarg, nhat);
    return opening_radians * tmpl_Rad_to_Deg;
}
