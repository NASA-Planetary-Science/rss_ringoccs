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
 *      1.) This function sets the tau->error_occured Boolean to true on      *
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

/*  Include all necessary headers.                                            */
#include <stdlib.h>
#include <math.h>

#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_special_functions.h>

#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  As a side note, the expression #var is used in the following macros. This *
 *  takes the input of the macro, which we're calling var for all three, and  *
 *  treats it as a string literal.                                            */

/*  Use this macro to save on repetitive code. It checks if tau->var is NULL, *
 *  attempts to malloc memory for tau->var if it is, and then checks to see   *
 *  if malloc failed.                                                         */
#define MALLOC_TAU_VAR(var)                                                    \
    /*  Check if the variable is not NULL. It should be at the start.        */\
    if (tau->var != NULL)                                                      \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_strdup(                                      \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"                         \
            "\r"#var" is not NULL. It is likely you've already set the data\n" \
            "\rfor this tau object. Returning.\n"                              \
        );                                                                     \
        return;                                                                \
    }                                                                          \
                                                                               \
    /*  Allocate memory for the variable.                                    */\
    tau->var = malloc(sizeof(*tau->var) * tau->arr_size);                      \
                                                                               \
    /*  Check if malloc failed.                                              */\
    if (tau->var == NULL)                                                      \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_strdup(                                      \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"                         \
            "\rMalloc failed and returned NULL for "#var". Returning.\n\n"     \
        );                                                                     \
        return;                                                                \
    }
/*  End of the __MALLOC_TAU_VAR__ macro.                                      */

/*  Use this macro to save on repetitive code. It is for checking that all of *
 *  of the values of a given member in the tau object are non-negative.       */
#define TAU_CHECK_NON_NEGATIVE(var)                                            \
    /*  Use tmpl_Double_Array_Min to compute the minimum value and check if  */\
    /*  it is negative. Return error if it is.                               */\
    if (tmpl_Double_Array_Min(tau->var, tau->arr_size) < 0.0)                  \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_strdup(                                      \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"                         \
            "\r"#var" has negative valued entries. Returning.\n"               \
        );                                                                     \
        return;                                                                \
    }
/*  End of the __TAU_CHECK_NON_NEGATIVE__ macro.                              */

/*  Use this macro to save on repetitive code. It is for checking that all of *
 *  of the values of a given member in the tau object fall within [-2pi, 2pi].*
 *  Note, it implicitly has min and max defined. These are declared at the    *
 *  top of the rssringoccs_Copy_DLP_Data_To_Tau function.                     */
#define TAU_CHECK_TWO_PI(var)                                                  \
    /*  Compute the minimum and maximum of var.                              */\
    tmpl_Double_Array_MinMax(tau->var, tau->arr_size, &min, &max);             \
                                                                               \
    /*  Check if var falls within the interval [-2pi, 2pi].                  */\
    if ((min < -tmpl_Two_Pi) || (max > tmpl_Two_Pi))                           \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_strdup(                                      \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"                         \
            "\r"#var" has values outside of [-2pi, 2pi]. Returning.\n"         \
        );                                                                     \
        return;                                                                \
    }
/*  End of the __TAU_CHECK_TWO_PI__ macro.                                    */

/*  Function for copying the relevant DLP data to a tau object.               */
void rssringoccs_Copy_DLP_Data_To_Tau(rssringoccs_DLPObj *dlp,
                                      rssringoccs_TAUObj *tau)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    unsigned long n;

    /*  These variables are used in the __TAU_CHECK_TWO_PI__ macro. Do not    *
     *  remove their declarations.                                            */
    double min, max;

    /*  If the tau pointer is NULL, we can't access it. Return.               */
    if (tau == NULL)
        return;

    /*  Exit the function if tau has its error_occurred member set to true.   */
    if (tau->error_occurred)
        return;

    /*  If the dlp pointer is NULL, raise and error and return.               */
    if (dlp == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "Input dlp pointer is NULL. Returning.\n"
        );
        return;
    }

    /*  If the dlp has its error_occurred member set to true, set the tau     *
     *  pointer's error_occurred member to true and exit.                     */
    if (dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "\rInput dlp pointer has the error_occurred member set to True.\n"
        );
        return;
    }

    /*  The arr_size for tau is the same as the dlp. Set this.                */
    tau->arr_size = dlp->arr_size;

    /*  If arr_size is less than 2, we can't do any processing and we can't   *
     *  compute dx_km. Return with error.                                     */
    if (tau->arr_size <= 1)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "\rInput arrays have less than 2 points. It is impossible to\n"
            "\rperform reconstrunction. Returning.\n"
        );
        return;
    }

    /*  The __MALLOC_TAU_VAR__ macro ends with an if statement and so has     *
     *  braces {}. Because of this, we do not need a semi-colon at the end.   *
     *  This macro allocates memory for the members of the tau object and     *
     *  checks for errors.                                                    */
    MALLOC_TAU_VAR(rho_km_vals)
    MALLOC_TAU_VAR(phi_rad_vals)
    MALLOC_TAU_VAR(f_sky_hz_vals)
    MALLOC_TAU_VAR(rho_dot_kms_vals)
    MALLOC_TAU_VAR(raw_tau_threshold_vals)
    MALLOC_TAU_VAR(B_rad_vals)
    MALLOC_TAU_VAR(D_km_vals)
    MALLOC_TAU_VAR(t_oet_spm_vals)
    MALLOC_TAU_VAR(t_ret_spm_vals)
    MALLOC_TAU_VAR(t_set_spm_vals)
    MALLOC_TAU_VAR(rho_corr_pole_km_vals)
    MALLOC_TAU_VAR(rho_corr_timing_km_vals)
    MALLOC_TAU_VAR(phi_rl_rad_vals)
    MALLOC_TAU_VAR(p_norm_vals)
    MALLOC_TAU_VAR(phase_rad_vals)
    MALLOC_TAU_VAR(rx_km_vals)
    MALLOC_TAU_VAR(ry_km_vals)
    MALLOC_TAU_VAR(rz_km_vals)

    /*  Loop through the entries of all of the pointers and set the nth value *
     *  of a tau member to the nth value of the corresponding dlp member.     */
    for (n=0; n<dlp->arr_size; ++n)
    {
        tau->rho_km_vals[n] = dlp->rho_km_vals[n];
        tau->phi_rad_vals[n] = dlp->phi_rad_vals[n];
        tau->B_rad_vals[n] = dlp->B_rad_vals[n];
        tau->D_km_vals[n] = dlp->D_km_vals[n];
        tau->f_sky_hz_vals[n] = dlp->f_sky_hz_vals[n];
        tau->rho_dot_kms_vals[n] = dlp->rho_dot_kms_vals[n];
        tau->t_oet_spm_vals[n] = dlp->t_oet_spm_vals[n];
        tau->t_ret_spm_vals[n] = dlp->t_ret_spm_vals[n];
        tau->t_set_spm_vals[n] = dlp->t_set_spm_vals[n];
        tau->rho_corr_pole_km_vals[n] = dlp->rho_corr_pole_km_vals[n];
        tau->rho_corr_timing_km_vals[n] = dlp->rho_corr_timing_km_vals[n];
        tau->phi_rl_rad_vals[n] = dlp->phi_rl_rad_vals[n];
        tau->p_norm_vals[n] = dlp->p_norm_vals[n];
        tau->raw_tau_threshold_vals[n] = dlp->raw_tau_threshold_vals[n];
        tau->rx_km_vals[n] = dlp->rx_km_vals[n];
        tau->ry_km_vals[n] = dlp->ry_km_vals[n];
        tau->rz_km_vals[n] = dlp->rz_km_vals[n];

        /*  The phase needs to be negated due to mathematical conventions.    */
        tau->phase_rad_vals[n] = -dlp->phase_rad_vals[n];
    }

    /*  Compute dx_km from the first and zeroth entries of rho_km_vals.       */
    tau->dx_km = tau->rho_km_vals[1] - tau->rho_km_vals[0];

    /*  Check if dx_km is a legal value.                                      */
    if (tau->dx_km == 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "\rdx_km = 0. rho_km_vals[1] - rho_km_vals[0] = 0. Returning.\n"
        );
        return;
    }

    /*  dx_km may be negative if this is an ingress occultation. To check if  *
     *  res is a legal value, compare it with twice the absolute value of     *
     *  dx_km. To avoid floating round-off error (which has happened the      *
     *  Cassini team, hence this edit) set the value to 1.99 instead of 2.0.  */
    else if (tau->res < 1.99 * fabs(tau->dx_km))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "\rResolution is less than twice the sample space.\n"
            "\rThis will result in an inaccurate reconstruction. Returning.\n"
        );
        return;
    }

    /*  The following members of tau should be non-negative for all entries.  *
     *  The __TAU_CHECK_NON_NEGATIVE__ contains an if statement and ends with *
     *  braces {}. Hence we do not need a semi-colon at the end.              */
    TAU_CHECK_NON_NEGATIVE(rho_km_vals)
    TAU_CHECK_NON_NEGATIVE(D_km_vals)
    TAU_CHECK_NON_NEGATIVE(f_sky_hz_vals)
    TAU_CHECK_NON_NEGATIVE(p_norm_vals)

    /*  Check that the following variables for angles fall within [-2pi, 2pi].*
     *  Like the other two macros, the __TAU_CHECK_TWO_PI__ macro ends with   *
     *  braces {} so we do not need a semi-colon at the end of these lines.   */
    TAU_CHECK_TWO_PI(B_rad_vals)
    TAU_CHECK_TWO_PI(phi_rad_vals)
    TAU_CHECK_TWO_PI(phase_rad_vals)
}
/*  End of rssringoccs_Copy_DLP_Data_To_Tau.                                  */
