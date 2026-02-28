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
 *                      rss_ringoccs_tau_check_resolution                     *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks for errors in the resolution found in a tau object.            *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Resolution                                      *
 *  Purpose:                                                                  *
 *      Checks for a few simple errors found in the resolution.               *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          A pointer to a tau object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Is_Inf:                                               *
 *              Checks if a double is +/- infinity.                           *
 *          tmpl_Double_Is_NaN:                                               *
 *              Checks if a double is NaN (Not-a-Number).                     *
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Notes:                                                                    *
 *      1.) If the error_occurred Boolean was previously set to true, this    *
 *          function does nothing and skips all checks.                       *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          A detailed description of the geometry of a ring occultation can  *
 *          be found here, and this includes the definition of resolution.    *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          tmpl_Double_Abs declared here.                                    *
 *  3.) rss_ringoccs_tau.h:                                                   *
 *          Tau definition and prototype for the function given here.         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       April 11, 2025                                                *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/02/28: Ryan Maguire                                                  *
 *      Fixed typos in doc-string, added status message when verbose is on.   *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Abs declared here.                                            */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Checks the resolution in a tau object for errors.                         */
void rssringoccs_Tau_Check_Resolution(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (tau->verbose)
        puts("\r\tTAU: Checking resolution for simple errors...");

    /*  The resolution should be a real number. Check for NaN (Not-a-Number). */
    if (tmpl_Double_Is_NaN(tau->resolution_km))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Resolution\n\n"
            "\rdx_km is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Resolution should also be finite. Treat infinity as an error.         */
    if (tmpl_Double_Is_Inf(tau->resolution_km))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Resolution\n\n"
            "\rdx_km is infinite.\n\n";

        return;
    }

    /*  A valid resolution must be positive.                                  */
    if (tau->resolution_km <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Resolution\n\n"
            "\rInput resolution is not positive.\n\n";

        return;
    }

    /*  dx_km may be negative if this is an ingress occultation. To check if  *
     *  the resolution is a legal value, compare it with twice the magnitude  *
     *  of dx_km. To avoid floating round-off error (which has happened to    *
     *  the Cassini team, hence this edit) set the value to 1.99 times dx.    */
    if (tau->resolution_km < 1.99 * tmpl_Double_Abs(tau->dlp->dx_km))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Resolution\n\n"
            "\rResolution is less than twice the sample spacing.\n"
            "\rThis will result in an inaccurate reconstruction.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Resolution.                                  */
