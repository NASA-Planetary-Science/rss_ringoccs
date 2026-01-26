/******************************************************************************
 *                                   LICENSE                                  *
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
 *                    rss_ringoccs_tau_check_keywords                         *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks all of the keywords specified by a user for the tau object.    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Keywords                                        *
 *  Purpose:                                                                  *
 *      Runs an error check on a rssringoccs_TAUObj pointer.                  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      src/tau/                                                              *
 *          rssringoccs_Tau_Check_Allan_Deviation:                            *
 *              Checks the Allan deviation (sigma) for errors.                *
 *          rssringoccs_Tau_Check_Resolution:                                 *
 *              Checks the requested resolution for simple errors.            *
 *          rssringoccs_Tau_Check_Eccentricity:                               *
 *              Checks the eccentricity for errors.                           *
 *          rssringoccs_Tau_Check_Periapse:                                   *
 *              Checks the angle of periapse for errors.                      *
 *          rssringoccs_Tau_Check_Range:                                      *
 *              Checks the processing range for very simple errors.           *
 *  Method:                                                                   *
 *      Check the keywords for the Tau object. There are several:             *
 *                                                                            *
 *          1.) Resolution.                                                   *
 *          2.) Allan deviation (sigma).                                      *
 *          3.) Periapse.                                                     *
 *          4.) Eccentricity.                                                 *
 *          5.) Processing range.                                             *
 *                                                                            *
 *      Each property is checked for basic errors (NaN, infinity, etc.).      *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL. Nothing is done in this case.      *
 *                                                                            *
 *      2.) If the error_occurred Boolean is true, nothing is done.           *
 *                                                                            *
 *      3.) This function sets the error_occurred Boolean to true on failure. *
 *          Inspect this (and the error_message) after calling this function. *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) rss_ringoccs_tau.h:                                                   *
 *          Header providing the function prototype and TAUObj type.          *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 4, 2021                                               *
 ******************************************************************************/

/*  Include the necessary header files.                                       */
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking the keyword arguments of a tau object.              */
void rssringoccs_Tau_Check_Keywords(rssringoccs_TAUObj * const tau)
{
    /*  If tau is NULL, then there is nothing to do. Return.                  */
    if (!tau)
        return;

    /*  Do not attempt to access data in a Tau object that previously had an  *
     *  error occur. Inspect the Boolean and abort if necessary.              */
    if (tau->error_occurred)
        return;

    rssringoccs_Tau_Check_Allan_Deviation(tau);
    rssringoccs_Tau_Check_Resolution(tau);
    rssringoccs_Tau_Check_Eccentricity(tau);
    rssringoccs_Tau_Check_Periapse(tau);
    rssringoccs_Tau_Check_Range(tau);
}
/*  End of rssringoccs_Tau_Check_Keywords.                                    */
