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
 *                       rss_ringoccs_dlp_check_geometry                      *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks for errors in the geometry data found in a dlp object.         *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Geometry                                        *
 *  Purpose:                                                                  *
 *      Runs all of the geometry checks on a DLP object.                      *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          A pointer to a dlp object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      src/dlp/                                                              *
 *          rssringoccs_DLP_Check_Azimuth_Angle:                              *
 *              Checks the azimuth angle "phi" for errors.                    *
 *          rssringoccs_DLP_Check_Displacement:                               *
 *              Checks the displacement "dx_km" for errors.                   *
 *          rssringoccs_DLP_Check_Opening_Angle:                              *
 *              Checks the opening angle "B" for errors.                      *
 *          rssringoccs_DLP_Check_Ring_Distance:                              *
 *              Checks the ring distance "D" for errors.                      *
 *          rssringoccs_DLP_Check_Ring_Radius:                                *
 *              Checks the ring radius "rho" for errors.                      *
 *  Notes:                                                                    *
 *      1.) If the error_occurred Boolean was previously set to true, this    *
 *          function does nothing and skips all checks.                       *
 *                                                                            *
 *      2.) This function checks for NULL pointers. Nothing is done if the    *
 *          input is NULL.                                                    *
 *                                                                            *
 *      3.) Each function that is called sets the "error_occurred" Boolean to *
 *          true on error. Inspect this after calling this routine.           *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          A detailed description of the geometry of a ring occultation can  *
 *          be found here, including the definitions of phi, B, D, and rho.   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) rss_ringoccs_dlp.h:                                                   *
 *          DLP object definition given here and helper routines provided.    *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       April 11, 2025                                                *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/02/27: Ryan Maguire                                                  *
 *      Made this a DLP error check, moved to the dlp folder.                 *
 ******************************************************************************/

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Function for checking DLP parameters for possible errors.                 */
void rssringoccs_DLP_Check_Geometry(rssringoccs_DLPObj * const dlp)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!dlp)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (dlp->error_occurred)
        return;

    /*  If the input is valid, run all of the geometry checks. Each function  *
     *  will check the error_occurred Boolean before proceeding and then set  *
     *  this Boolean to true on error. Because of this we do not need to      *
     *  repeatedly check the error_occurred Boolean between function calls,   *
     *  instead we may simply chain all of the function calls together.       */
    rssringoccs_DLP_Check_Azimuth_Angle(dlp);
    rssringoccs_DLP_Check_Displacement(dlp);
    rssringoccs_DLP_Check_Opening_Angle(dlp);
    rssringoccs_DLP_Check_Ring_Distance(dlp);
    rssringoccs_DLP_Check_Ring_Radius(dlp);
}
/*  End of rssringoccs_DLP_Check_Geometry.                                    */
