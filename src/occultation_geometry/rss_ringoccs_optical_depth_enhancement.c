#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>

double
rssringoccs_Double_Optical_Depth_Enhancement(const double b, const double phi)
{
    const double b_eff = rssringoccs_Double_Effective_Ring_Opening(b, phi);
    return 1.0 / tmpl_Double_Tan(b_eff);
}
