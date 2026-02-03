/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Typedef for the Tau object and functions working with the struct.     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       July 21, 2023                                                 *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_TAU_H
#define RSS_RINGOCCS_TAU_H

/*  TMPL_RESTRICT given here.                                                 */
#include <libtmpl/include/tmpl_config.h>

/*  Booleans (True / False) provided here.                                    */
#include <libtmpl/include/tmpl_bool.h>

/*  Complex numbers and functions provided here.                              */
#include <libtmpl/include/types/tmpl_complex_double.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  DLP object is typedef'd here.                                             */
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>

/*  Tau object typedef provided here.                                         */
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>

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
    const double * TMPL_RESTRICT const,       /*  The array r - r0.           */
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

extern const rssringoccs_FresnelNewtonTransform
rssringoccs_newton_transform_table[7];

extern const rssringoccs_FresnelTransform
rssringoccs_newton_interp_transform_table[6];

extern rssringoccs_FresnelNewtonTransform
rssringoccs_Tau_Select_Newton_Transform(rssringoccs_TAUObj * const tau);

extern rssringoccs_FresnelTransform
rssringoccs_Tau_Select_Newton_Interp_Transform(rssringoccs_TAUObj * const tau);

extern void
rssringoccs_Tau_Check_Allan_Deviation(rssringoccs_TAUObj * const tau);

extern void rssringoccs_Tau_Check_Eccentricity(rssringoccs_TAUObj * const tau);
extern void rssringoccs_Tau_Check_Periapse(rssringoccs_TAUObj * const tau);
extern void rssringoccs_Tau_Check_Range(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Azimuth_Angle                                   *
 *  Purpose:                                                                  *
 *      Checks the phi_deg_vals array in a tau object for common errors.      *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Azimuth_Angle(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Core_Data                                       *
 *  Purpose:                                                                  *
 *      Checks the core pointers in a tau objects for errors.                 *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The tau object we are checking.                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Core_Data(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Data_Range                                      *
 *  Purpose:                                                                  *
 *      Checks if there is enough data to process a data set.                 *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The tau object we are checking.                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Data_Range(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Displacement                                    *
 *  Purpose:                                                                  *
 *      Checks the dx_km value in a tau object for common errors.             *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Displacement(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Geometry                                        *
 *  Purpose:                                                                  *
 *      Checks a Tau object for possible errors.                              *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Geometry(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Opening_Angle                                   *
 *  Purpose:                                                                  *
 *      Checks the B_deg_vals array in a tau object for common errors.        *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Opening_Angle(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Resolution                                      *
 *  Purpose:                                                                  *
 *      Checks the resolution_km value in a tau object for common errors.     *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Resolution(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Ring_Distance                                   *
 *  Purpose:                                                                  *
 *      Checks the D_km_vals array in a tau object for common errors.         *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Ring_Distance(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Check_Ring_Radius                                     *
 *  Purpose:                                                                  *
 *      Checks the rho_km_vals array in a tau object for common errors.       *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Check_Ring_Radius(rssringoccs_TAUObj * const tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Create_From_DLP                                       *
 *  Purpose:                                                                  *
 *      Creates a Tau object using the data in a DLP objects.                 *
 *  Arguments:                                                                *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *      res (double):                                                         *
 *          The requested resolution, in kilometers.                          *
 *  Outputs:                                                                  *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The output Tau object.                                            *
 *  Notes:                                                                    *
 *      malloc is called to allocate memory for the output pointer and the    *
 *      members in the Tau object. If malloc fails to allocate memory for the *
 *      tau pointer, then NULL is returned. If malloc fails to allocate       *
 *      memory for any of the members of the Tau objects, the error_occurred  *
 *      Boolean is set to true. Check both of these before accessing data.    *
 ******************************************************************************/
extern rssringoccs_TAUObj *
rssringoccs_Tau_Create_From_DLP(const rssringoccs_DLPObj *dlp, double res);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Copy_DLP_Data                                         *
 *  Purpose:                                                                  *
 *      Copies data and computes variables for a Tau object from DLP data.    *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object. DLP members will be copied here.                  *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Copy_DLP_Data(
    const rssringoccs_DLPObj * TMPL_RESTRICT const dlp,
    rssringoccs_TAUObj * TMPL_RESTRICT const tau
);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Copy_DLP_Members                                      *
 *  Purpose:                                                                  *
 *      Copies members that dlp and tau objects have in common from a dlp     *
 *      object in to a tau object.                                            *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object. DLP members will be copied here.                  *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Copy_DLP_Members(rssringoccs_TAUObj *tau,
                                 const rssringoccs_DLPObj *dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Compute_Data_From_DLP_Members                         *
 *  Purpose:                                                                  *
 *      Computes several Tau variables from the given DLP data.               *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object.                                                   *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object.                                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Compute_Data_From_DLP_Members(rssringoccs_TAUObj *tau,
                                              const rssringoccs_DLPObj *dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Destroy                                               *
 *  Purpose:                                                                  *
 *      Frees all data associated with a Tau object.                          *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj **):                                          *
 *          The Tau object that is to be destroyed.                           *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Destroy(rssringoccs_TAUObj **tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Destroy_Members                                       *
 *  Purpose:                                                                  *
 *      Frees all data associated with a Tau object without destroying the    *
 *      actual Tau pointer.                                                   *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to be freed.                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Destroy_Members(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Init                                                  *
 *  Purpose:                                                                  *
 *      Initialize a tau struct so that it's members are NULL.                *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to initialized.                  *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      Functions that allocate memory for a Tau object check if the members  *
 *      are NULL before working with them. If the pointer is not NULL it is   *
 *      assumed memory was successfully allocated for the variable. Always    *
 *      call this function on a new tau object first.                         *
 ******************************************************************************/
extern void rssringoccs_Tau_Init(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Malloc_Members                                        *
 *  Purpose:                                                                  *
 *      Allocates memory for the Tau variables.                               *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to be allocated memory.          *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      It is assumed tau->arr_size has been set and that tau is not NULL.    *
 ******************************************************************************/
extern void rssringoccs_Tau_Malloc_Members(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Default_Values                                    *
 *  Purpose:                                                                  *
 *      Sets the default values for several parameters.                       *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose values are to be set.                        *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Set_Default_Values(rssringoccs_TAUObj* tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Psi_Type                                          *
 *  Purpose:                                                                  *
 *      Sets the reconstruction method corresponding to a string.             *
 *  Arguments:                                                                *
 *      psitype (const char *):                                               *
 *          A string containing the requested reconstruction method.          *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose psinum is to be set.                         *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Psi_Type(const char *psitype, rssringoccs_TAUObj* tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Range_From_String                                 *
 *  Purpose:                                                                  *
 *      Sets the reconstruction range corresponding to a string.              *
 *  Arguments:                                                                *
 *      range (const char *):                                                 *
 *          A string containing the requested range of reconstruction.        *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose range is to be set.                          *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Range_From_String(const char *range,
                                      rssringoccs_TAUObj* tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Window_Type                                       *
 *  Purpose:                                                                  *
 *      Sets the window function in a Tau object corresponding to a string.   *
 *  Arguments:                                                                *
 *      wtype (const char *):                                                 *
 *          A string containing the requested window funtion.                 *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose window function is to be set.                *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Window_Type(const char *wtype, rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Keywords(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Occ_Type(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Reset_Window(const rssringoccs_TAUObj * TMPL_RESTRICT const tau,
                             double * TMPL_RESTRICT const x_arr,
                             double * TMPL_RESTRICT const w_func,
                             const size_t n_pts,
                             const size_t center);

extern void
rssringoccs_Tau_Compute_Window(rssringoccs_TAUObj * TMPL_RESTRICT const tau,
                               double * TMPL_RESTRICT const w_func,
                               const size_t nw_pts,
                               const size_t center);

extern int
rssringoccs_Tau_Resize_Half_Window(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    double ** TMPL_RESTRICT const x_ptr,
    double ** TMPL_RESTRICT const w_ptr,
    const double width,
    const double two_dx,
    const size_t center
);

#endif
/*  End of include guard.                                                     */
