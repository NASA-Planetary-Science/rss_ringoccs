
#ifndef RSS_RINGOCCS_OCCULTATION_GEOMETRY_H
#define RSS_RINGOCCS_OCCULTATION_GEOMETRY_H

#include <libtmpl/include/tmpl.h>

extern double rssringoccs_Double_Light_Distance(double t0, double t1);

extern double
rssringoccs_Double_Ring_Opening_Angle(const double time,
                                      const tmpl_ThreeVectorDouble * const nhat,
                                      const char *const spacecraft,
                                      const char *const dsn,
                                      const char *const kernels,
                                      const char *const ref);

extern double
rssringoccs_Double_Optical_Depth_Enhancement(const double b, const double phi);

extern double
rssringoccs_Double_Effective_Ring_Opening(const double b, const double phi);
#endif
