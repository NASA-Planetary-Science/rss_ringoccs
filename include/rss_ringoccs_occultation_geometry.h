

#ifndef RSS_RINGOCCS_OCCULTATION_GEOMETRY_H
#define RSS_RINGOCCS_OCCULTATION_GEOMETRY_H

#include <libtmpl/include/tmpl_euclidean_spatial_geometry.h>

extern void
rssringoccs_Calc_B_Deg(double *et_vals, char *spacecraft, char *dsn,
                       tmpl_ThreeVector nhat_p, char *kernels, char *ref);

#endif

