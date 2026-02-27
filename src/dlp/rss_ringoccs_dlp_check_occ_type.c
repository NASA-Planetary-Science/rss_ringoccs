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
 *                     rss_ringoccs_dlp_check_occ_type                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Check the data stored in an rssringoccs_DLPObj pointer and determines *
 *      if it is an ingress or egress occultation, or a chord.                *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Occ_Type:                                       *
 *  Purpose:                                                                  *
 *      Determination what type of occultation the DLP data contains.         *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          A pointer to a rssringoccs_DLPObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      1.) This function sets the dlp->error_occured Boolean to true on      *
 *          error. It is the user's responsibility to check that this Boolean *
 *          is false after using this function. Trying to access the pointers *
 *          in a rssringoccs_DLPObj may result in a segmentation fault        *
 *          otherwise.                                                        *
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
 *  5.) rss_ringoccs_reconstruction.h:                                        *
 *          The rssringoccs_DLPObj is defined here and the function           *
 *          prototypes for reconstruction are found here as well.             *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 2, 2021                                               *
 ******************************************************************************/

#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>
#include <stddef.h>

/*  Function for determining the type of occultation of a dlp object.         */
void rssringoccs_DLP_Check_Occ_Type(rssringoccs_DLPObj * const dlp)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double min;
    double max;
    size_t n;

    /*  Check if the dlp pointer is NULL, returning if it is.                 */
    if (dlp == NULL)
        return;

    /*  If the dlp object has its error_occurred member set to true, do not   *
     *  do any computations and return.                                       */
    if (dlp->error_occurred)
        return;

    /*  Check that the pointers we need to access are not NULL. If they are,  *
     *  the user forgot to copy the relevant data from the DLP object or      *
     *  prematurely destroyed/free'd the data from dlp.                       */
    if (dlp->rho_dot_kms_vals == NULL)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\rInput rho_dot_kms_vals is NULL.\n\n";

        return;
    }

    /*  Is arr_size is less than 2 we can't do any processing. Return error.  */
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
            "\r\tdrho/dt has positive and negative values.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction on each part.\n\n";

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
            "\r\tdrho/dt has zero valued elements.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction on each part.\n\n";

        return;
    }

    /*  If dx_km is negative and rho_dot_kms_vals is positive, there is most  *
     *  likely an error. Rather than assuming what the occultation is and     *
     *  proceeding with fingers cross, return an error. The user should fix   *
     *  DLP data so the dx_km is positive.                                    */
    if ((dlp->dx_km < 0.0) && (min > 0.0))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Occ_Type\n\n"
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\trho_km_vals is decreasing yet rho_dot_kms_vals\n"
            "\r\tis positiive. Check DLP data for errors.\n\n";

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
     *  rho_dot_kms_vals. rssringoccs_Reverse_Double_Array is found in        *
     *  in rss_ringoccs_math.h.                                               */
    else if (dlp->dx_km < 0.0)
    {
        tmpl_Double_Array_Reverse(dlp->rho_km_vals,      dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->phi_deg_vals,     dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->B_deg_vals,       dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->D_km_vals,        dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->rho_dot_kms_vals, dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->t_oet_spm_vals,   dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->t_ret_spm_vals,   dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->t_set_spm_vals,   dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->phi_rl_deg_vals,  dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->rx_km_vals,       dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->ry_km_vals,       dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->rz_km_vals,       dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->rho_corr_pole_km_vals, dlp->arr_size);
        tmpl_Double_Array_Reverse(dlp->rho_corr_timing_km_vals, dlp->arr_size);

        for(n = 0; n < dlp->arr_size; ++n)
            dlp->rho_dot_kms_vals[n] =
                tmpl_Double_Abs(dlp->rho_dot_kms_vals[n]);

        dlp->dx_km = tmpl_Double_Abs(dlp->dx_km);
    }
}
/*  End of rssringoccs_DLP_Check_Occ_Type.                                    */
