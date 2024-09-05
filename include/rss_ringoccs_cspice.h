

#ifndef RSS_RINGOCCS_CSPICE_H
#define RSS_RINGOCCS_CSPICE_H

extern void kclear_c(void);
extern void furnsh_c(const char *kernels);
extern void
spkpos_c(const char *targ, double et, const char *ref, const char *abcorr,
         const char *obs, double ptarg[3], double *lt);

#endif
