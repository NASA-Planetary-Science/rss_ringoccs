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



extern void
rssringoccs_Tau_Check_Data(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Keywords(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Occ_Type(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Reset_Window(double * TMPL_RESTRICT const x_arr,
                             double * TMPL_RESTRICT const w_func,
                             double dx,
                             double width,
                             size_t nw_pts,
                             rssringoccs_WindowFunction fw);

extern int
rssringoccs_Tau_Resize_Half_Window(rssringoccs_TAUObj * const tau,
                                   double ** const x_ptr,
                                   double ** const w_ptr,
                                   double width,
                                   double two_dx,
                                   size_t center);

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
