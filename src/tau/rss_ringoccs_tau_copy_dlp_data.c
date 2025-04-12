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
 *                   rss_ringoccs_copy_dlp_data_to_tau                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provides the code needed for copying all relevant data from a         *
 *      rssringoccs_DLPObj pointer into a rssringoccs_TAUObj pointer.         *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Copy_DLP_Data_To_Tau:                                     *
 *  Purpose:                                                                  *
 *      Copies the relevant data from a rssringoccs_DLPObj pointer to a       *
 *      rssringoccs_TAUObj pointer.                                           *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          A pointer to a rssringoccs_DLPObj.                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Method:                                                                   *
 *      Malloc memory for the members of the rssringoccs_TAUObj pointer and   *
 *      loop through the indices copying the values from the                  *
 *      rssringoccs_DLPObj pointer.                                           *
 *  NOTES:                                                                    *
 *      1.) This function sets the tau->error_occurred Boolean to true on     *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_TAUObj may result in a segmentation fault        *
 *          otherwise.                                                        *
 *      2.) The tau object does NOT steal the reference to the pointers in    *
 *          the dlp pointer. Destroying the tau object does NOT free the      *
 *          memory in the dlp object, and vice-versa. To avoid memory leaks   *
 *          one must destroy both the dlp and tau objects when done with them.*
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) stdlib.h:                                                             *
 *          C standard library header. Used for the NULL macro and malloc.    *
 *  2.) rss_ringoccs_bool.h:                                                  *
 *          Header file containing rssringoccs_Bool, and True and False.      *
 *  3.) rss_ringoccs_math.h:                                                  *
 *          Header file containing various math routines. This header file    *
 *          provides compatibility between the C89/C90 and C99 math.h files.  *
 *          The C99 version is a superset of the C89 one. This header file    *
 *          aliases various functions if C99 is available, and defines the    *
 *          missing ones otherwise.                                           *
 *  4.) rss_ringoccs_string.h:                                                *
 *          Header file containing routines for manipulating strings. The     *
 *          rssringoccs_strdup function is defined here. strdup is a function *
 *          that comes with POSIX but is not part of the C standard. Because  *
 *          of this, rss_ringoccs provides an implementation of this that     *
 *          only uses C89/C90 compliant code.                                 *
 *  5.) rss_ringoccs_calibration.h:                                           *
 *          Contains the typedef for the rssringoccs_DLPObj structure.        *
 *  6.) rss_ringoccs_reconstruction.h:                                        *
 *          The rssringoccs_TAUObj is defined here and the function           *
 *          prototypes for reconstruction are found here as well.             *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 1, 2021                                               *
 ******************************************************************************/

/*  Include all necessary headers.                                            */
#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Function for copying the relevant DLP data to a tau object.               */
void
rssringoccs_Tau_Copy_DLP_Data(
    const rssringoccs_DLPObj * TMPL_RESTRICT const dlp,
    rssringoccs_TAUObj * TMPL_RESTRICT const tau
)
{
    /*  If the tau pointer is NULL, we can't access it. Return.               */
    if (!tau)
        return;

    /*  Exit the function if tau has its error_occurred member set to true.   */
    if (tau->error_occurred)
        return;

    /*  If the dlp pointer is NULL, raise and error and return.               */
    if (!dlp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Copy_DLP_Data\n\n"
            "\rInput dlp pointer is NULL.\n\n";

        return;
    }

    /*  If the dlp has its error_occurred member set to true, set the tau     *
     *  pointer's error_occurred member to true and exit.                     */
    if (dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Copy_DLP_Data\n\n"
            "\rInput dlp error_occurred set to True.\n";

        return;
    }

    /*  If arr_size is less than 2, we can't do any processing and we can't   *
     *  compute dx_km. Return with error.                                     */
    if (dlp->arr_size <= 1)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Copy_DLP_Data\n\n"
            "\rdlp->arr_size < 2, arrays have less than 2 points.\n\n";

        return;
    }

    /*  The arr_size for tau is the same as the dlp. Set this.                */
    tau->arr_size = dlp->arr_size;

    /*  Allocate memory for the tau variables.                                */
    rssringoccs_Tau_Malloc_Members(tau);

    /*  Check for errors before continuing. The previous function calls       *
     *  malloc and it may possibly fail to allocate memory.                   */
    if (tau->error_occurred)
        return;

    /*  Several variables for Tau are exactly the same as the DLP ones. These *
     *  can be copied verbatim into the newly allocated memory for Tau.       */
    rssringoccs_Tau_Copy_DLP_Members(tau, dlp);

    /*  Other variables need to be computed from the DLP data. Do this.       */
    rssringoccs_Tau_Compute_Data_From_DLP_Members(tau, dlp);

    /*  Lastly, compute dx from the first and zeroth entries of rho_km_vals.  */
    tau->dx_km = tau->rho_km_vals[1] - tau->rho_km_vals[0];
}
/*  End of rssringoccs_Tau_Copy_DLP_Data.                                     */
