/*  Include guard for this file to prevent including this twice.              */
#ifndef _RSS_RINGOCCS_PPMPLOT_H_
#define _RSS_RINGOCCS_PPMPLOT_H_

#include <stdio.h>

extern void
rssringoccs_Color(unsigned char red, unsigned char green,
                  unsigned char blue, FILE *fp);

extern void rssringoccs_RGB_Scaled_Gradient(double val, FILE *fp);

extern void
rssringoccs_RGB_Linear_Gradient(double val, double min, double max, FILE *fp);

#endif
