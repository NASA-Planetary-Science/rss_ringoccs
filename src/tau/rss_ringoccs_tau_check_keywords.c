/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                    rss_ringoccs_tau_check_keywords                         *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks all of the keywords specified by a user for the tau object.    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Keywords:                                       *
 *  Purpose:                                                                  *
 *      Runs an error check on a rssringoccs_TAUObj pointer.                  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Method:                                                                   *
 *      Check that the keywords have legal values, setting errors if not.     *
 *  NOTES:                                                                    *
 *      1.) This function sets the tau->error_occured Boolean to true on      *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_TAUObj may result in a segmentation fault        *
 *          otherwise.                                                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) stdlib.h:                                                             *
 *          C standard library header. Used for the NULL macro.               *
 *  2.) rss_ringoccs_bool.h:                                                  *
 *          Header file containing rssringoccs_Bool, and True and False.      *
 *  3.) rss_ringoccs_string.h:                                                *
 *          Header file containing routines for manipulating strings. The     *
 *          rssringoccs_strdup function is defined here. strdup is a function *
 *          that comes with POSIX but is not part of the C standard. Because  *
 *          of this, rss_ringoccs provides an implementation of this that     *
 *          only uses C89/C90 compliant code.                                 *
 *  4.) rss_ringoccs_reconstruction.h:                                        *
 *          The rssringoccs_TAUObj is defined here and the function           *
 *          prototypes for reconstruction are found here as well.             *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in config_librssringoccs.sh uses gcc and has the  *
 *  -pedantic and -std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 4, 2021                                               *
 ******************************************************************************/

/*  Include the necessary header files.                                       */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking the keyword arguments of a tau object.              */
void rssringoccs_Tau_Check_Keywords(rssringoccs_TAUObj *tau)
{
    /*  If tau is NULL we can't access its members. Return.                   */
    if (!tau)
        return;

    /*  Check if the error_occurred member was set to true before proceeding. */
    if (tau->error_occurred)
        return;

    rssringoccs_Tau_Check_Resolution(tau);

    /*  The Allen deviation must be positive.                                 */
    if (tau->sigma <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput sigma (Allen deviation) is not positive.\n\n";

        return;
    }

    /*  The periapse is allowed to be between -2pi and 2pi, inclusive.        */
    if (tau->periapse < -tmpl_Double_Two_Pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse less than -2pi.\n\n";

        return;
    }

    if (tau->periapse > tmpl_Double_Two_Pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse greater than 2pi.\n\n";

        return;
    }

    /*  Eccentricity can be positive or zero. No negative values allowed.     */
    if (tau->eccentricity < 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput eccentricity is negative.\n\n";

        return;
    }

    /*  Lastly, check the requested range values.                             */
    if (tau->rng_list[0] < 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rStarting value for range is negative.\n\n";

        return;
    }

    if (tau->rng_list[0] > tau->rng_list[1])
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rStarting value for range is greater than final value.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Keywords.                                    */
