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
 *                     rss_ringoccs_tau_check_data_range                      *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Check the data stored in an rssringoccs_TAUObj pointer and determines *
 *      if there is enough data to perform diffraction correction.            *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Data_Range                                      *
 *  Purpose:                                                                  *
 *      Runs an error check on a rssringoccs_TAUObj pointer.                  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_string.h:                                                        *
 *          tmpl_String_Duplicate:                                            *
 *              Copies a string and stores it in another string.              *
 *  Method:                                                                   *
 *      Loop through the requested reconstruction region and see if there is  *
 *      enough data to the left and right of each point to perform the        *
 *      inversion. Given a point rho, there should be rho-W/2 to rho+W/2 of   *
 *      data where W is the window width stored in tau->w_km_vals.            *
 *  Notes:                                                                    *
 *      1.) This function sets the tau->error_occurred Boolean to true on     *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_TAUObj may result in a segmentation fault        *
 *          otherwise.                                                        *
 *      2.) No data is freed if an error occurs. The caller must do this.     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans (True and False).                  *
 *  2.) compat/tmpl_cast.h:                                                   *
 *          Macros for casting with compatibility for both C and C++.         *
 *  4.) tmpl_string.h:                                                        *
 *          Provides POSIX functions like strdup in a portable manner.        *
 *  6.) rss_ringoccs_tau.h:                                                   *
 *          Header file where the rssringoccs_TAUObj type is provided.        *
 *  7.) stddef.h:                                                             *
 *          Standard library header file containing size_t typedef.           *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 1, 2021                                               *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2024/12/23: Ryan Maguire                                                  *
 *      Cleaned up code, added some comments, fixed includes.                 *
 ******************************************************************************/

/*  Booleans (True / False) provided here.                                    */
#include <libtmpl/include/tmpl_bool.h>

/*  Macros for C vs. C++ compatibility with casting data types.               */
#include <libtmpl/include/compat/tmpl_cast.h>

/*  Portable version of strdup provided here.                                 */
#include <libtmpl/include/tmpl_string.h>

/*  rssringoccs_TAUObj typedef given here, and the function prototype.        */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef found here.                                                */
#include <stddef.h>

/*  Function for checking the data range of a rssringoccs_TAUObj pointer.     */
void rssringoccs_Tau_Check_Data_Range(rssringoccs_TAUObj * const tau)
{
    /*  C89 does not allow mixed code. All declarations are made at the top.  *
     *  Declare a variable to check what the maximum index is.                */
    size_t current_max;

    /*  Variables for indexing the for loop.                                  */
    size_t n, start, end;

    /*  Variable for the reciprocal of 2.0 * dx.                              */
    double rcpr_two_dx;

    /*  And a variable for the number of points in a window.                  */
    size_t nw_pts;

    /*  Check that the input pointer is not NULL before trying to access it.  */
    if (!tau)
        return;

    /*  If an error occurred before calling this function, abort.             */
    if (tau->error_occurred)
        return;

    /*  Check that tau->dx_km is not zero to avoid divide-by-zero.            */
    if (tau->dx_km <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Data_Range\n\n"
            "\rtau->dx_km is not positive. Returning.\n"
        );

        return;
    }

    /*  Set the values of the user-requested min and max indices.             */
    start = tau->start;
    end = start + tau->n_used;

    /*  Set the rcpr_two_dx value from the tau object. Division is more       *
     *  expensive computationally than multiplication, so we store the        *
     *  reciprocal of 2 * dx and compute with that.                           */
    rcpr_two_dx = 0.5 / tau->dx_km;

    /*  Loop through every point, check window width, and ensure you have     *
     *  enough data to the left and right for data processing.                */
    for (n = start; n < end; ++n)
    {
        /*  Compute the number of points needed in a window.                  */
        nw_pts = TMPL_CAST(tau->w_km_vals[n] * rcpr_two_dx, size_t);

        /*  If n - nw_pts is negative, then the window goes beyond the        *
         *  available data. Since our integer variables are declared as       *
         *  size_t, which is an unsigned data type, the expression n - nw_pts *
         *  can't be negative. To avoid error, we simply check if n < nw_pts. *
         *  If it is, we raise an error and return.                           */
        if (n < nw_pts)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_String_Duplicate(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_Data_Range\n\n"
                "\rNot enough data to perform diffraction correction. The\n"
                "\rrequested region has points with a window width that go\n"
                "\rbeyond the minimum radius available. Returning.\n"
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
            tau->error_message = tmpl_String_Duplicate(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_Data_Range\n\n"
                "\rNot enough data to perform diffraction correction. The\n"
                "\rrequested region has points with a window width that go\n"
                "\rbeyond the maximum radius available. Returning.\n"
            );

            return;
        }
    }
    /*  End of for loop computing the number of points for the windows.       */
}
/*  End of rssringoccs_Tau_Check_Data_Range.                                  */
