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
 *                     rss_ringoccs_tau_check_occ_type                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Check the data stored in an rssringoccs_TAUObj pointer and determines *
 *      if it is an ingress or egress occultation, or a chord.                *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Check_Occ_Type:                                       *
 *  Purpose:                                                                  *
 *      Determination what type of occultation the tau data contains.         *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          A pointer to a rssringoccs_TAUObj.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
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
 *          The rssringoccs_TAUObj is defined here and the function           *
 *          prototypes for reconstruction are found here as well.             *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 2, 2021                                               *
 ******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Function for determining the type of occultation of a tau object.         */
void rssringoccs_Tau_Check_Occ_Type(rssringoccs_TAUObj *tau)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double min;
    double max;
    size_t n;

    /*  Check if the tau pointer is NULL, returning if it is.                 */
    if (tau == NULL)
        return;

    /*  If the tau object has its error_occurred member set to true, do not   *
     *  do any computations and return.                                       */
    if (tau->error_occurred)
        return;

    /*  Check that the pointers we need to access are not NULL. If they are,  *
     *  the user forgot to copy the relevant data from the DLP object or      *
     *  prematurely destroyed/free'd the data from tau.                       */
    if (tau->rho_dot_kms_vals == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\rInput rho_dot_kms_vals is NULL. Returning.\n"
        );
        return;
    }

    /*  Is arr_size is less than 2 we can't do any processing. Return error.  */
    if (tau->arr_size <= 1)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\rInput arrays have less than 2 points. It is impossible to\n"
            "\rperform reconstrunction. Returning.\n"
        );
        return;
    }

    /*  If dx_km is zero, return with error. We can't check the occ type.     */
    if (tau->dx_km == 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\rdx_km is zero. Impossible to determine occ type. Returning.\n"
        );
        return;
    }

    /*  Compute the minimum and maximum of rho_dot_kms_vals.                  */
    tmpl_Double_Array_MinMax(tau->rho_dot_kms_vals, tau->arr_size, &min, &max);

    /*  If rho_dot_kms_vals has both negative and positive values, then the   *
     *  occultation is a chord-occ and we can't continue, or there is an      *
     *  error in the DLP data. In either case, return with error.             */
    if ((min < 0.0) && (max > 0.0))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\r\tdrho/dt has positive and negative values.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return;
    }

    /*  If there are entries of rho_dot_kms_vals that are zero, again it is   *
     *  likely that the data comes from a chord occultation and the file was  *
     *  improperly split into egress and ingress portions. Return error.      */
    else if ((min == 0.0) || (max == 0.0))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\r\tdrho/dt has zero valued elements.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return;
    }

    /*  If rho_dot_kms_vals is negative but dx_km is positive, then we have   *
     *  an ingress occultation and need to compute with the absolute value    *
     *  of rho_dot_kms_vals. Compute this and store in the tau object.        */
    else if ((tau->dx_km > 0.0) && (max < 0.0))
    {
        for(n=0; n<tau->arr_size; ++n)
            tau->rho_dot_kms_vals[n] = fabs(tau->rho_dot_kms_vals[n]);
    }

    /*  If dx_km is negative and rho_dot_kms_vals is positive, there is most  *
     *  likely an error. Rather than assuming what the occultation is and     *
     *  proceeding with fingers cross, return an error. The user should fix   *
     *  DLP data so the dx_km is positive.                                    */
    else if ((tau->dx_km < 0.0) && (min > 0.0))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\trho_km_vals is decreasing yet rho_dot_kms_vals\n"
            "\r\tis positiive. Check DLP data for errors.\n\n"
        );
        return;
    }

    /*  If dx_km is negative, and if rho_dot_kms_vals is not zero or positive *
     *  as ruled out by the previous if-else-then statements, then we can     *
     *  safely assume ingress with the data decreasing in radius. Reverse the *
     *  data to be increasing in radius and compute the absolute value of     *
     *  rho_dot_kms_vals. rssringoccs_Reverse_Double_Array is found in        *
     *  in rss_ringoccs_math.h.                                               */
    else if (tau->dx_km < 0.0)
    {
        tmpl_Double_Array_Reverse(tau->rho_km_vals,      tau->arr_size);
        tmpl_Double_Array_Reverse(tau->phi_deg_vals,     tau->arr_size);
        tmpl_Double_Array_Reverse(tau->B_deg_vals,       tau->arr_size);
        tmpl_Double_Array_Reverse(tau->D_km_vals,        tau->arr_size);
        tmpl_Double_Array_Reverse(tau->rho_dot_kms_vals, tau->arr_size);
        tmpl_Double_Array_Reverse(tau->t_oet_spm_vals,   tau->arr_size);
        tmpl_Double_Array_Reverse(tau->t_ret_spm_vals,   tau->arr_size);
        tmpl_Double_Array_Reverse(tau->t_set_spm_vals,   tau->arr_size);
        tmpl_Double_Array_Reverse(tau->phi_rl_deg_vals,  tau->arr_size);
        tmpl_Double_Array_Reverse(tau->rx_km_vals,       tau->arr_size);
        tmpl_Double_Array_Reverse(tau->ry_km_vals,       tau->arr_size);
        tmpl_Double_Array_Reverse(tau->rz_km_vals,       tau->arr_size);
        tmpl_Double_Array_Reverse(tau->rho_corr_pole_km_vals, tau->arr_size);
        tmpl_Double_Array_Reverse(tau->rho_corr_timing_km_vals, tau->arr_size);

        for(n=0; n<tau->arr_size; ++n)
            tau->rho_dot_kms_vals[n] = tmpl_Double_Abs(tau->rho_dot_kms_vals[n]);

        tau->dx_km = -tau->dx_km;
    }
}
/*  End of rssringoccs_Tau_Check_Occ_Type.                                    */
