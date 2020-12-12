/*  Include guard for this file to prevent including this twice.              */
#ifndef _RSS_RINGOCCS_PPMPLOT_H_
#define _RSS_RINGOCCS_PPMPLOT_H_

#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

extern void
rssringoccs_Color(unsigned char red, unsigned char green,
                  unsigned char blue, FILE *fp);

extern void rssringoccs_RGB_Scaled_Gradient(double val, FILE *fp);

extern void
rssringoccs_RGB_Linear_Gradient(double val, double min, double max, FILE *fp);

extern void
rssringoccs_Easy_Complex_Plots(
    const char *func_name,
    rssringoccs_ComplexDouble(*f)(rssringoccs_ComplexDouble),
    unsigned int x_size, unsigned int y_size,
    const double x_min, const double x_max,
    const double y_min, const double y_max
);

extern void
rssringoccs_Easy_Real_Plots(const char *func_name, double (*f)(double),
                            unsigned int x_size, unsigned int y_size,
                            const double x_min, const double x_max,
                            const double y_min, const double y_max);

#endif
