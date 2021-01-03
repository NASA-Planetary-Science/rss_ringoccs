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
 *                      rss_ringoccs_check_tau_data                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks the relevant points in a rssringoccs_TAUObj pointer to see if  *
 *      they are NULL. Sets an error if they are.                             *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Check_Tau_Data:                                           *
 *  Purpose:                                                                  *
 *      Runs an error check on a rssringoccs_TAUObj pointer.                  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Method:                                                                   *
 *      Check the relevant members to see if they are NULL.                   *
 *  NOTES:                                                                    *
 *      1.) This function sets the tau->error_occured Boolean to true on      *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_TAUObj may result in a segmentation fault        *
 *          otherwise.                                                        *
 *      2.) This function is called by all of the diffraction correction      *
 *          routines at the start as a safety net for avoiding segmentation   *
 *          faults from occuring and crashing the program.                    *
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
 *  Date:       January 1, 2021                                               *
 ******************************************************************************/

/*  Include the necessary header files.                                       */
#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Function for checking the pointers in a rssringoccs_TAUObj pointer.       */
void rssringoccs_Check_Tau_Data(rssringoccs_TAUObj *tau)
{
    /*  First, check that the actual rssringoccs_TAUObj pointer is not NULL.  */
    if (tau == NULL)
        return;

    /*  If the tau->error_occurred member is true, return.                    */
    if (tau->error_occurred)
        return;

    /*  If any of these pointers are NULL, set the error_occurred member to   *
     *  true, set an error message, and return.                               */
    if (tau->T_in == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Check_Tau_Data\n\n"
            "\rInput tau has T_in set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->T_out == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Check_Tau_Data\n\n"
            "\rInput tau has T_out set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->rho_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has rho_km_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->F_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has F_km_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->phi_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has phi_rad_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->kd_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has kd_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->B_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has B_rad_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->D_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has D_km_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->w_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has w_km_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    /*  If of if-statements checking if pointers are NULL.                    */
}
/*  End of rssringoccs_Check_Tau_Data.                                        */
