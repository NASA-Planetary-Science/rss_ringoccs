
#ifndef RSS_RINGOCCS_OCCULTATION_GEOMETRY_H
#define RSS_RINGOCCS_OCCULTATION_GEOMETRY_H

#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/types/tmpl_vec3_double.h>

extern double
rssringoccs_Double_Ring_Opening_Angle_Deg(
    const double time,
    const tmpl_ThreeVectorDouble * TMPL_RESTRICT const nhat,
    const char * TMPL_RESTRICT const spacecraft,
    const char * TMPL_RESTRICT const dsn,
    const char * TMPL_RESTRICT const ref
);

extern double
rssringoccs_Double_Elevation_Angle_Deg(
    const double time,
    const char * TMPL_RESTRICT const target,
    const char * TMPL_RESTRICT const observer,
    const char * TMPL_RESTRICT const ref
);

#endif
