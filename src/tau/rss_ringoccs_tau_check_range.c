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
 *                        rss_ringoccs_tau_check_range                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks the range member of a Tau object for simple errors.            *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Range                                           *
 *  Purpose:                                                                  *
 *      Checks the range member of a Tau object for simple errors.            *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          A pointer to the Tau object that we are checking.                 *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Is_Inf:                                               *
 *              Checks if a double is infinity.                               *
 *          tmpl_Double_Is_NaN:                                               *
 *              Checks if a double is Not-a-Number.                           *
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Method:                                                                   *
 *      Check the range for simple errors:                                    *
 *          1.) Check that range[0] and range[1] are not NaN.                 *
 *          2.) Make sure that range[0] and range[1] are not infinity.        *
 *          3.) range[0] and range[1] should be positive.                     *
 *          4.) The range should be ordered, range[0] < range[1].             *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers and it checks if the       *
 *          error_occurred Boolean is set to true. Nothing is done in either  *
 *          of these cases.                                                   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          Header providing Is_NaN and Is_Inf functions.                     *
 *  3.) rss_ringoccs_tauobj.h:                                                *
 *          Tau object typedef provided here.                                 *
 *  4.) stdio.h:                                                              *
 *          Standard library header file providing the puts function.         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 23, 2026                                              *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/02/27: Ryan Maguire                                                  *
 *      Simple clean-up, added doc-string.                                    *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Is_NaN and tmpl_Double_Is_Inf declared here.                  */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Function for checking if the input range is valid.                        */
void rssringoccs_Tau_Check_Range(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (tau->dlp->verbose)
        puts("\r\tTAU: Checking range array for simple errors...");

    /*  Range values should be real numbers. Check for NaN (Not-a-Number).    */
    if (tmpl_Double_Is_NaN(tau->range[0]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[0] is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Same check for the next value in the range list.                      */
    if (tmpl_Double_Is_NaN(tau->range[1]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[1] is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  The range should also be finite. Treat infinity as an error.          */
    if (tmpl_Double_Is_Inf(tau->range[0]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[0] is infinite.\n\n";

        return;
    }

    /*  Same check for the next value in the range list.                      */
    if (tmpl_Double_Is_Inf(tau->range[1]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[1] is infinite.\n\n";

        return;
    }

    /*  The lower range should be positive.                                   */
    if (tau->range[0] < 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rStarting value for range is negative.\n\n";

        return;
    }

    /*  Lastly, the range list should be in increasing order.                 */
    if (tau->range[0] > tau->range[1])
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rStarting value for range is greater than final value.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Range.                                       */
