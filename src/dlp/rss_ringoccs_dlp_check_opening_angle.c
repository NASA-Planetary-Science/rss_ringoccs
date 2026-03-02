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
 *                    rss_ringoccs_dlp_check_opening_angle                    *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks for errors in the ring opening angle found in a dlp object.    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Opening_Angle                                   *
 *  Purpose:                                                                  *
 *      Checks for a few common errors found in the ring opening angle of     *
 *      a dlp object.                                                         *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          A pointer to a dlp object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Array_MinMax:                                         *
 *              Computes the min and max of a double array.                   *
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Notes:                                                                    *
 *      1.) It is assumed that the ring opening angle has been allocated      *
 *          memory and the data has been initialized. If B_deg_vals is NULL,  *
 *          the error_occurred Boolean will be set to True.                   *
 *      2.) This function checks for NULL pointers. If dlp is NULL, nothing   *
 *          is done. If B_deg_vals is NULL, an error message is set.          *
 *      3.) If the error_occurred Boolean was previously set to true, this    *
 *          function does nothing and skips all checks.                       *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          A detailed description of the geometry of a ring occultation can  *
 *          be found here, and this includes the definition of B.             *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          tmpl_Double_Array_MinMax declared here.                           *
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

/*  tmpl_Double_Array_MinMax declared here, computes min and max of an array. */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Checks the ring opening angle in a dlp object for simple errors.          */
void rssringoccs_DLP_Check_Opening_Angle(rssringoccs_DLPObj * const dlp)
{
    /*  Variables for the min and max of the B_deg_vals array.                */
    double min, max;

    /*  If the input is NULL there is nothing to be done.                     */
    if (!dlp)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (dlp->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (dlp->verbose)
        puts("\r\tDLP: Checking opening angle for simple errors...");

    /*  This function should only be called after the B_deg_vals array was    *
     *  allocated memory and the data initialized. NULL pointers are hence    *
     *  treated as errors.                                                    */
    if (!dlp->B_deg_vals)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Opening_Angle\n\n"
            "\rB_deg_vals is NULL.\n\n";

        return;
    }

    /*  Compute the minimum and maximum of B_deg_vals.                        */
    tmpl_Double_Array_MinMax(dlp->B_deg_vals, dlp->arr_size, &min, &max);

    /*  The data should be reduced mod 360 since the angle is in degrees. If  *
     *  min(B_deg_vals) < -360, then an error likely occurred.                */
    if (min < -360.0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Opening_Angle\n\n"
            "\rOpening angle (B) has values less than -360 degrees.\n\n";

        return;
    }

    /*  Similarly, the maximum angle should not be larger than 360.           */
    if (max > 360.0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Opening_Angle\n\n"
            "\rOpening angle (B) has angles greater than +360 degrees.\n\n";
    }
}
/*  End of rssringoccs_DLP_Check_Opening_Angle.                               */
