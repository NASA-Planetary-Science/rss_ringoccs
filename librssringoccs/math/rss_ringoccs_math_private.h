#ifndef __RSS_RINGOCCS_MATH_PRIVATE_H__
#define __RSS_RINGOCCS_MATH_PRIVATE_H__

/*  Auxilliary functions for computing sine and cosine.                       */
extern float rssringoccs_do_sinf(float x);
extern double rssringoccs_do_sin(double x);
extern long double rssringoccs_do_sinl(long double x);

extern float rssringoccs_do_cosf(float x);
extern double rssringoccs_do_cos(double x);
extern long double rssringoccs_do_cosl(long double x);

extern float rssringoccs_sinf_table(unsigned int n);
extern double rssringoccs_sin_table(unsigned int n);
extern long double rssringoccs_sinl_table(unsigned int n);

extern float rssringoccs_cosf_table(unsigned int n);
extern double rssringoccs_cos_table(unsigned int n);
extern long double rssringoccs_cosl_table(unsigned int n);

#endif
