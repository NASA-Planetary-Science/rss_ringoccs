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
 *                      rss_ringoccs_tau_check_periapse                       *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks the periapse (theta) for simple errors.                        *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Periapse                                        *
 *  Purpose:                                                                  *
 *      Checks the periapse in a Tau object for simple errors.                *
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
 *      Check the periapse for simple errors:                                 *
 *          1.) Check that sigma is not NaN.                                  *
 *          2.) Make sure that sigma is not infinity.                         *
 *          3.) The periapse should be between -2 pi and 2 pi.                *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers and it checks if the       *
 *          error_occurred Boolean is set to true. Nothing is done in either  *
 *          of these cases.                                                   *
 *  References:                                                               *
 *      1.) Gresh, Marouf, Tyler, Rosen, Simpson (July 1988)                  *
 *          Voyager Radio Occultation by Uranus' Rings                        *
 *          Icarus Vol. 78, Pages 131-168.                                    *
 *                                                                            *
 *          This paper covers Fresnel optics with elliptical rings. The       *
 *          elliptic code in rss_ringoccs is based on this.                   *
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

/*  Header file with the Tau definition.                                      */
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_Tau_Check_Periapse(rssringoccs_TAUObj * const tau);

/*  Function for checking if the input periapse is valid.                     */
void rssringoccs_Tau_Check_Periapse(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (tau->verbose)
        puts("\r\tTAU: Checking periapse for simple errors...");

    /*  The periapse should be a real number. Check for NaN (Not-a-Number).   */
    if (tmpl_Double_Is_NaN(tau->periapse))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Periapse\n\n"
            "\rperiapse is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Periapse should also be finite. Treat infinity as an error.           */
    if (tmpl_Double_Is_Inf(tau->periapse))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Periapse\n\n"
            "\rperiapse is infinite.\n\n";

        return;
    }

    /*  The periapse is allowed to be between -2pi and 2pi, inclusive.        */
    if (tau->periapse < -tmpl_double_two_pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse less than -2 pi.\n\n";

        return;
    }

    /*  Same check, peripase should be bounded by +2 pi.                      */
    if (tau->periapse > tmpl_double_two_pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse greater than +2 pi.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Periapse.                                    */
