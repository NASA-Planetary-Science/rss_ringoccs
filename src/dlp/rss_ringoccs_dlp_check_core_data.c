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
 *                      rss_ringoccs_tau_check_core_data                      *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks the core pointers in a tau object to see if they are NULL.     *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Core_Data                                       *
 *  Purpose:                                                                  *
 *      Runs an error check on a tau object, ensuring the core arrays are not *
 *      not NULL.                                                             *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          A pointer to a tau object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      Check the relevant members to see if they are NULL.                   *
 *  Notes:                                                                    *
 *      1.) This function sets the tau->error_occurred Boolean to true on     *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a tau object may result in a segmentation fault otherwise.     *
 *      2.) No data is freed if an error occurs. The caller must do this.     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans (True and False).                  *
 *  2.) rss_ringoccs_tau.h:                                                   *
 *          Header file where the rssringoccs_TAUObj type is provided.        *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 1, 2021                                               *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2025/04/15: Ryan Maguire                                                  *
 *      More clean up, removed tmpl_string.h dependency, moved to tau folder. *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Macro for checking the data in a tau object. This is to save repetitive   *
 *  code, it simply checks if a certain pointer in tau is NULL. The #var      *
 *  preprocessor directive treats var as a string literal. Note that since    *
 *  this macro ends with braces, we do not need semi-colons when calling it.  */
#define RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(var)                                 \
    if (!tau->var)                                                             \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message =                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Check_Core_Data\n\n"                          \
            "\rInput tau has "#var" set to NULL.\n\n";                         \
        return;                                                                \
    }
/*  End of RSSRINGOCCS_TAU_CHECK_DATA_MEMBER macro.                           */

/*  Function for checking the core pointers in a tau object.                  */
void rssringoccs_Tau_Check_Core_Data(rssringoccs_TAUObj * const tau)
{
    /*  Check that the input pointer is not NULL before trying to access it.  */
    if (!tau)
        return;

    /*  If an error occurred before calling this function, abort.             */
    if (tau->error_occurred)
        return;

    /*  The following are the core pointers in a tau object, the data that is *
     *  explicitly needed for diffraction correction. If any of them are NULL *
     *  this is to be treated as an error. Note that the macro defined above, *
     *  RSSRINGOCCS_TAU_CHECK_DATA_MEMBER, contains an if-then with braces.   *
     *  Because of this we do not need semi-colons at the end of these lines. */
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(T_in)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(T_out)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(rho_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(F_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(phi_deg_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(k_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(rho_dot_kms_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(B_deg_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(D_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(w_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(rx_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(ry_km_vals)
    RSSRINGOCCS_TAU_CHECK_DATA_MEMBER(rz_km_vals)
}
/*  End of rssringoccs_Tau_Check_Core_Data.                                   */

/*  Undefine everything in case someone wants to #include this file.          */
#undef RSSRINGOCCS_TAU_CHECK_DATA_MEMBER
