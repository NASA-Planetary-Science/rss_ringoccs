/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RECONSTRUCTION_H
#define RSS_RINGOCCS_RECONSTRUCTION_H

#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef is given here.                                             */
#include <stdlib.h>

/*  Fresnel transform. Inputs are Tau data, window data, the number of points *
 *  in the window, and the center of the window (its index in the data).      */
typedef void
(*rssringoccs_FresnelNewtonTransform)(
    rssringoccs_TAUObj * TMPL_RESTRICT const, /*  Reconstruction data.        */
    const double * TMPL_RESTRICT const,       /*  Tapering / window array.    */
    size_t,                                   /*  Number of points in window. */
    size_t                                    /*  Index for center of window. */
);

typedef void (*rssringoccs_FresnelLegendreTransform)(
    rssringoccs_TAUObj * TMPL_RESTRICT const, /*  Reconstruction data.        */
    const double * TMPL_RESTRICT const,       /*  The array (r - r0) / D.     */
    const double * TMPL_RESTRICT const,       /*  Tapering / window array.    */
    const double * TMPL_RESTRICT const,       /*  Polynomial coefficients.    */
    size_t,                                   /*  Number of points in window. */
    size_t                                    /*  Index for center of window. */
);

typedef void (*rssringoccs_FresnelTransform)(
    rssringoccs_TAUObj * TMPL_RESTRICT const, /*  Reconstruction data.        */
    const double * TMPL_RESTRICT const,       /*  The array (r - r0) / D.     */
    const double * TMPL_RESTRICT const,       /*  Tapering / window array.    */
    size_t,                                   /*  Number of points in window. */
    size_t                                    /*  Index for center of window. */
);

extern void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau);

/*  Functions that compute the Fresnel Transform on a TAUObj instance.        */
extern void
rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Legendre(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_SimpleFFT(rssringoccs_TAUObj *tau);

#endif
