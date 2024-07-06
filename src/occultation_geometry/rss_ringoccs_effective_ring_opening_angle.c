#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_occultation_geometry.h>

double
rssringoccs_Double_Effective_Ring_Opening(const double b, const double phi)
{
    const double b_rad = tmpl_Deg_to_Rad * b;
    return tmpl_Double_Arctan2(tmpl_Double_Tan(b_rad), tmpl_Double_Cosd(phi));
}
