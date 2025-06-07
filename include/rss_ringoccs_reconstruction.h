/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RECONSTRUCTION_H
#define RSS_RINGOCCS_RECONSTRUCTION_H

#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef is given here.                                             */
#include <stddef.h>

extern void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau);

/*  Functions that compute the Fresnel Transform on a TAUObj instance.        */
extern void
rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Legendre(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj * const tau);

extern void
rssringoccs_Diffraction_Correction_Polynomial_Newton(
    rssringoccs_TAUObj * const tau
);

extern void
rssringoccs_Diffraction_Correction_SimpleFFT(rssringoccs_TAUObj *tau);

extern void rssringoccs_Diffraction_Correction(rssringoccs_TAUObj *tau);

#endif
