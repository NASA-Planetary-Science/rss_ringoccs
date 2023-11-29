/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RECONSTRUCTION_H
#define RSS_RINGOCCS_RECONSTRUCTION_H

/*  Various functions, complex variables, and more found here.                */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef is given here.                                             */
#include <stdlib.h>

/*  Fresnel transform. Inputs are Tau data, window data, the number of points *
 *  in the window, and the center of the window (its index in the data).      */
typedef void
(*rssringoccs_FresT)(rssringoccs_TAUObj *, const double *, size_t, size_t);

extern void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Data_Range(rssringoccs_TAUObj *dlp);

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
rssringoccs_Tau_Reset_Window(double *x_arr, double *w_func, double dx,
                             double width, size_t nw_pts,
                             rssringoccs_Window_Function fw);

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
