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
 *                      rss_ringoccs_dlp_check_occ_type                       *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Check the data stored in an rssringoccs_DLPObj pointer and determines *
 *      if it is an ingress or egress occultation, or a chord.                *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Occ_Type                                        *
 *  Purpose:                                                                  *
 *      Determination what type of occultation the DLP data contains.         *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          A pointer to a rssringoccs_DLPObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      src/dlp/                                                              *
 *          rssringoccs_DLP_Reverse_Occultation:                              *
 *              Swaps the order of the arrays in a DLP object.                *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Array_MinMax:                                         *
 *              Find the min and max of an array.                             *
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Notes:                                                                    *
 *      1.) This function sets the dlp->error_occurred Boolean to true on     *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_DLPObj may result in a segmentation fault        *
 *          otherwise.                                                        *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          Header providing math routines.                                   *
 *  3.) rss_ringoccs_dlp.h:                                                   *
 *          DLP object typedef provided here.                                 *
 *  4.) stdio.h:                                                              *
 *          Standard library header file providing the puts function.         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 2, 2021                                               *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2025/04/15: Ryan Maguire                                                  *
 *      Moved routine from the reconstruction folder into the tau directory.  *
 *  2026/02/27: Ryan Maguire                                                  *
 *      Made this a DLP error check, moved to the dlp folder.                 *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Array_MinMax declared here, computes min and max of an array. */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Function for determining the type of occultation of a dlp object.         */
void rssringoccs_DLP_Check_Occ_Type(rssringoccs_DLPObj * const dlp)
{
    /*  Variables for the min and max of the rho_dot_kms_vals array.          */
    double min, max;

    /*  Variable for indexing the rho_dot_kms_vals array.                     */
    size_t n;

    /*  Check if the dlp pointer is NULL, returning if it is.                 */
    if (!dlp)
        return;

    /*  If the dlp object has its error_occurred member set to true, do not   *
     *  do any computations and return.                                       */
    if (dlp->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (dlp->verbose)
        puts("\r\tDLP: Checking occultation type (ingress vs. egress)...");

    /*  Check that the pointers we need to access are not NULL. If they are,  *
     *  the user forgot to copy the relevant data for the DLP object or       *
     *  prematurely destroyed / free'd the data from dlp.                     */
    if (!dlp->rho_dot_kms_vals)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rInput rho_dot_kms_vals is NULL.\n\n";

        return;
    }

    /*  If arr_size is less than 2 we can't do any processing. Return error.  */
    if (dlp->arr_size <= 1)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rInput arrays have less than 2 points.\n\n";

        return;
    }

    /*  If dx_km is zero, return with error. We can't check the occ type.     */
    if (dlp->dx_km == 0.0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rdx_km is zero. Impossible to determine occultation type.\n\n";

        return;
    }

    /*  Compute the minimum and maximum of rho_dot_kms_vals.                  */
    tmpl_Double_Array_MinMax(dlp->rho_dot_kms_vals, dlp->arr_size, &min, &max);

    /*  If rho_dot_kms_vals has both negative and positive values, then the   *
     *  occultation is a chord-occ and we can't continue, or there is an      *
     *  error in the DLP data. In either case, return with error.             */
    if ((min < 0.0) && (max > 0.0))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rdrho/dt has positive and negative values.\n"
            "\rYour input file is probably a chord occultation.\n"
            "\rDiffraction Correction can only be performed for\n"
            "\rone event at a time. That is, ingress or egress.\n\n"
            "\rTO CORRECT THIS:\n"
            "\r\tSplit the input into two parts: Ingress and Engress\n"
            "\r\tand perform diffraction correction on each part.\n\n";

        return;
    }

    /*  If there are entries of rho_dot_kms_vals that are zero, again it is   *
     *  likely that the data comes from a chord occultation and the file was  *
     *  improperly split into egress and ingress portions. Return error.      */
    if ((min == 0.0) || (max == 0.0))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rdrho/dt has zero valued elements.\n"
            "\rYour input file is probably a chord occultation.\n"
            "\rDiffraction Correction can only be performed for\n"
            "\rone event at a time. That is, ingress or egress.\n\n"
            "\rTO CORRECT THIS:\n"
            "\r\tSplit the input into two parts: Ingress and Engress\n"
            "\r\tand perform diffraction correction on each part.\n\n";

        return;
    }

    /*  If dx_km is negative and rho_dot_kms_vals is positive, there is most  *
     *  likely an error. Rather than assuming what the occultation is and     *
     *  proceeding with fingers crossed, return an error. The user should fix *
     *  the DLP data so that dx_km is positive.                               */
    if ((dlp->dx_km < 0.0) && (min > 0.0))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rrho_km_vals is decreasing yet rho_dot_kms_vals\n"
            "\ris positive. Check DLP data for errors.\n\n";

        return;
    }

    /*  If rho_dot_kms_vals is negative but dx_km is positive, then we have   *
     *  an ingress occultation and need to compute with the absolute value    *
     *  of rho_dot_kms_vals. Compute this and store in the dlp object.        */
    if ((dlp->dx_km > 0.0) && (max < 0.0))
    {
        for(n = 0; n < dlp->arr_size; ++n)
            dlp->rho_dot_kms_vals[n] =
                tmpl_Double_Abs(dlp->rho_dot_kms_vals[n]);
    }

    /*  If dx_km is negative, and if rho_dot_kms_vals is not zero or positive *
     *  as ruled out by the previous if-else-then statements, then we can     *
     *  safely assume ingress with the data decreasing in radius. Reverse the *
     *  data to be increasing in radius and compute the absolute value of     *
     *  rho_dot_kms_vals.                                                     */
    else if (dlp->dx_km < 0.0)
        rssringoccs_DLP_Reverse_Occultation(dlp);
}
/*  End of rssringoccs_DLP_Check_Occ_Type.                                    */
