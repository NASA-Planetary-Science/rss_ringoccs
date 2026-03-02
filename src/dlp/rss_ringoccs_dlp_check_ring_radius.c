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
 *                 rss_ringoccs_dlp_check_ring_intercept_point                *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks for errors in the ring radius found in a dlp object.           *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Ring_Radius                                     *
 *  Purpose:                                                                  *
 *      Checks for a few common errors found in the ring radius.              *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          A pointer to a dlp object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Array_Min:                                            *
 *              Computes the minimum of a double array.                       *
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Notes:                                                                    *
 *      1.) It is assumed that the ring radius has been allocated memory      *
 *          and the data has been initialized. If rho_km_vals is NULL, the    *
 *          error_occurred Boolean will be set to True.                       *
 *      2.) This function checks for NULL pointers. If dlp is NULL, nothing   *
 *          is done. If rho_km_vals is NULL, an error message is set.         *
 *      3.) If the error_occurred Boolean was previously set to true, this    *
 *          function does nothing and skips all checks.                       *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          A detailed description of the geometry of a ring occultation can  *
 *          be found here, and this includes the definition of rho.           *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          tmpl_Double_Array_Min declared here.                              *
 *  3.) rss_ringoccs_dlp.h:                                                   *
 *          DLP definition and prototype for the function given here.         *
 *  4.) stdio.h:                                                              *
 *          Standard library header file providing the puts function.         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       April 11, 2025                                                *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Array_Min declared here, computes the minimum of an array.    */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Checks the ring radius in a dlp object for errors.                        */
void rssringoccs_DLP_Check_Ring_Radius(rssringoccs_DLPObj * const dlp)
{
    /*  Variable for the minimum of the rho_km_vals array.                    */
    double min;

    /*  If the input is NULL there is nothing to be done.                     */
    if (!dlp)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (dlp->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (dlp->verbose)
        puts("\r\tDLP: Checking ring radius (rho) for simple errors...");

    /*  This function should only be called after the rho_km_vals array was   *
     *  allocated memory and the data initialized. NULL pointers are hence    *
     *  treated as errors.                                                    */
    if (!dlp->rho_km_vals)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Ring_Radius\n\n"
            "\rrho_km_vals is NULL.\n\n";

        return;
    }

    /*  Compute the minimum of rho_km_vals. This must be positive.            */
    min = tmpl_Double_Array_Min(dlp->rho_km_vals, dlp->arr_size);

    /*  rho_km_vals is a distance and it must be positive. Check for this.    */
    if (min <= 0.0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Ring_Radius\n\n"
            "\rRing radius (rho_km_vals) has non-positive values.\n\n";

        return;
    }
}
/*  End of rssringoccs_DLP_Check_Ring_Radius.                                 */
