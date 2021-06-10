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
 *                   rss_ringoccs_tau_check_data_range                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Check the data stored in an rssringoccs_TAUObj pointer and determines *
 *      if there is enough data to perform diffraction correction.            *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Check_Tau_Data_Range:                                     *
 *  Purpose:                                                                  *
 *      Runs an error check on a rssringoccs_TAUObj pointer.                  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Method:                                                                   *
 *      Loop through the requested reconstruction region and see if there is  *
 *      enough data to the left and right of each point to perform the        *
 *      inversion. Given a point rho, there should be rho-W/2 to rho+W/2 of   *
 *      data where W is the window width stored in tau->w_km_vals.            *
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
 *  Date:       January 1, 2021                                               *
 ******************************************************************************/

/*  Include the necessary header files.                                       */
#include <stdlib.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Function for checking the data range of a rssringoccs_TAUObj pointer.     */
void rssringoccs_Tau_Check_Data_Range(rssringoccs_TAUObj *tau)
{
    /*  As a side not, C89/C90 does not allow mixed code. All declaration     *
     *  must be made at the top of a code block.                              */

    /*  Declare a variable to check what the maximum index is.                */
    unsigned long current_max;

    /*  Variables for indexing the for loop.                                  */
    unsigned long n, start, end;

    /*  Variable for the reciprocal of 2.0 * dx.                              */
    double rcpr_two_dx;

    /*  And a variable for the number of points in a window.                  */
    unsigned long nw_pts;

    /*  Check that the input rssringoccs_TAUObj pointer is not NULL before    *
     *  trying to access it's members.                                        */
    if (tau == NULL)
        return;

    /*  Check that the rssringoccs_TAUObj pointer does not have its           *
     *  error_occurred member set to true.                                    */
    if (tau->error_occurred)
        return;

    /*  Check that tau->dx_km is not zero to avoid divide-by-zero.            */
    if (tau->dx_km <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Check_Tau_Data_Range\n\n"
            "\rtau->dx_km is not positive. Returning.\n"
        );
        return;
    }

    /*  Set the values of the user-requested min and max indices.             */
    start = tau->start;
    end   = start + tau->n_used;

    /*  Set the rcpr_two_dx value from the tau object. Division is more       *
     *  expensive computationally than multiplication, so we store the        *
     *  reciprocal of 2 * dx and compute with that.                           */
    rcpr_two_dx = 0.5 / tau->dx_km;

    /*  Loop through every point, check window width, and ensure you have     *
     *  enough data to the left and right for data processing.                */
    for (n = start; n < end; ++n)
    {
        /*  Compute the number of points needed in a window.                  */
        nw_pts = ((unsigned long)(tau->w_km_vals[n] * rcpr_two_dx));

        /*  If n - nw_pts is negative, then the window goes beyond the        *
         *  available data. Since our integer variables are declared as       *
         *  unsigned, n - nw_pts can't be negative. To avoid error, we simply *
         *  check if n < nw_pts. If it is, we raise an error and return.      */
        if (n < nw_pts)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Check_Tau_Data_Range\n\n"
                "\rNot enough data to perform diffraction correction. The\n"
                "\rrequested region has points with a window width that go\n"
                "\rbeyond the minimum radius you have. Returning.\n"
            );
            return;
        }

        /*  The largest radius needed for the window corresponds to the       *
         *  current point plus the number of points in the window.            */
        current_max = n + nw_pts;

        /*  If max_requested goes between the size of the array, we have      *
         *  illegal values. Return with error.                                */
        if (current_max > tau->arr_size)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Check_Tau_Data_Range\n\n"
                "\rNot enough data to perform diffraction correction. The\n"
                "\rrequested region has points with a window width that go\n"
                "\rbeyond the maximum radius you have. Returning.\n"
            );
            return;
        }
    }
    /*  End of for loop computing the number of points for the windows.       */
}
/*  End of rssringoccs_Check_Tau_Data_Range.                                  */
