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
 *      Typedef for the DLP object and functions working with the struct.     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       February 27, 2026                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_DLP_H
#define RSS_RINGOCCS_DLP_H

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  DLP object is typedef'd here.                                             */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Azimuth_Angle                                   *
 *  Purpose:                                                                  *
 *      Checks the phi_deg_vals array in a dlp object for common errors.      *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Azimuth_Angle(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Core_Data                                       *
 *  Purpose:                                                                  *
 *      Checks the core pointers in a dlp objects for errors.                 *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The dlp object we are checking.                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Core_Data(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Displacement                                    *
 *  Purpose:                                                                  *
 *      Checks the dx_km value in a dlp object for common errors.             *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Displacement(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Geometry                                        *
 *  Purpose:                                                                  *
 *      Checks a DLP object for possible errors in the geometry.              *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Geometry(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Opening_Angle                                   *
 *  Purpose:                                                                  *
 *      Checks the B_deg_vals array in a dlp object for common errors.        *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Opening_Angle(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Occ_Type                                        *
 *  Purpose:                                                                  *
 *      Checks the data and determines the occultation type.                  *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object whose values are to be checked.                    *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Occ_Type(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Ring_Distance                                   *
 *  Purpose:                                                                  *
 *      Checks the D_km_vals array in a dlp object for common errors.         *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Ring_Distance(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Check_Ring_Radius                                     *
 *  Purpose:                                                                  *
 *      Checks the rho_km_vals array in a dlp object for common errors.       *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          The DLP object to be checked.                                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Check_Ring_Radius(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Destroy                                               *
 *  Purpose:                                                                  *
 *      Frees all data associated with a DLP object.                          *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj **):                                          *
 *          The DLP object that is to be destroyed.                           *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Destroy(rssringoccs_DLPObj ** const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Destroy_Members                                       *
 *  Purpose:                                                                  *
 *      Frees all data associated with a DLP object without destroying the    *
 *      actual DLP pointer.                                                   *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object whose members are to be freed.                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Destroy_Members(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Init                                                  *
 *  Purpose:                                                                  *
 *      Initialize a dlp struct so that it's members are NULL.                *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object whose members are to initialized.                  *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      Functions that allocate memory for a DLP object check if the members  *
 *      are NULL before working with them. If the pointer is not NULL it is   *
 *      assumed memory was successfully allocated for the variable. Always    *
 *      call this function on a new dlp object first.                         *
 ******************************************************************************/
extern void rssringoccs_DLP_Init(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Malloc_Members                                        *
 *  Purpose:                                                                  *
 *      Allocates memory for the DLP variables.                               *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object whose members are to be allocated memory.          *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      It is assumed dlp->arr_size has been set and that dlp is not NULL.    *
 ******************************************************************************/
extern void rssringoccs_DLP_Malloc_Members(rssringoccs_DLPObj * const dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_DLP_Reverse_Occultation                                   *
 *  Purpose:                                                                  *
 *      Reverses the order of the arrays in a DLP object.                     *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          The DLP object whose members are to be flipped.                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_DLP_Reverse_Occultation(rssringoccs_DLPObj * const dlp);

extern rssringoccs_DLPObj *
rssringoccs_DLP_New_Reference(rssringoccs_DLPObj * const dlp);

extern void rssringoccs_DLP_Release(rssringoccs_DLPObj *dlp);

#endif
/*  End of include guard.                                                     */
