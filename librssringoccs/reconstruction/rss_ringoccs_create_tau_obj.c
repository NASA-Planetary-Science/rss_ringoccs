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
 *                      rss_ringoccs_create_tau_obj                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for creating an rssringoccs_TAUObj pointer.  *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Create_TAUObj:                                            *
 *  Purpose:                                                                  *
 *      Creates an rssringoccs_TAUObj pointer.                                *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj *):                                           *
 *          A pointer to a rssringoccs_DLPObj.                                *
 *      res (double):                                                         *
 *          The resolution for the diffraction correction, in kilometers.     *
 *  Output:                                                                   *
 *      tau (rssringoccs_TAUObj *).                                           *
 *          A pointer to a rssringoccs_TAUObj with the data copied over from  *
 *          the input dlp. The various default keywords are also set, and can *
 *          be changed later if desired.                                      *
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
 *  3.) rss_ringoccs_string.h:                                                *
 *          Header file containing routines for manipulating strings. The     *
 *          rssringoccs_strdup function is defined here. strdup is a function *
 *          that comes with POSIX but is not part of the C standard. Because  *
 *          of this, rss_ringoccs provides an implementation of this that     *
 *          only uses C89/C90 compliant code.                                 *
 *  4.) rss_ringoccs_calibration.h:                                           *
 *          Contains the typedef for the rssringoccs_DLPObj structure.        *
 *  5.) rss_ringoccs_reconstruction.h:                                        *
 *          The rssringoccs_TAUObj is defined here and the function           *
 *          prototypes for reconstruction are found here as well.             *
 *  6.) rss_ringoccs_special_functions.h:                                     *
 *          Header file containing definitions for various special functions. *
 *          Needed for the KBMD20NormEQ macro.                                *
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

/*  Include all necessary header files.                                       */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <stdlib.h>

/*  Function for allocating memory for a Tau object and setting the default   *
 *  values for all of the keywords.                                           */
rssringoccs_TAUObj *
rssringoccs_Create_TAUObj(rssringoccs_DLPObj *dlp, double res)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_TAUObj *tau;

    /*  Try to allocate memory for the tau pointer.                           */
    tau = malloc(sizeof(*tau));

    /*  Check if malloc failed.                                               */
    if (tau == NULL)
        return tau;

    /*  Before anything else, initialize all pointers of tau to NULL. When    *
     *  errors occur or when the tau object is destroyed, the pointers in tau *
     *  will attempt to be free'd. If they weren't malloc'd this will crash   *
     *  the program. The error routines and tau deallocator checks if these   *
     *  pointers are NULL before trying to free them. All pointers that are   *
     *  not NULL are assumed to be malloc'd pointers that can be free'd.      */
    tau->T_in = NULL;
    tau->T_out = NULL;
    tau->rho_km_vals = NULL;
    tau->F_km_vals = NULL;
    tau->phi_rad_vals = NULL;
    tau->k_vals = NULL;
    tau->f_sky_hz_vals = NULL;
    tau->rho_dot_kms_vals = NULL;
    tau->raw_tau_threshold_vals = NULL;
    tau->B_rad_vals = NULL;
    tau->D_km_vals = NULL;
    tau->w_km_vals = NULL;
    tau->t_oet_spm_vals = NULL;
    tau->t_ret_spm_vals = NULL;
    tau->t_set_spm_vals = NULL;
    tau->rho_corr_pole_km_vals = NULL;
    tau->rho_corr_timing_km_vals = NULL;
    tau->tau_threshold_vals = NULL;
    tau->phi_rl_rad_vals = NULL;
    tau->p_norm_vals = NULL;
    tau->phase_rad_vals = NULL;
    tau->wtype = NULL;
    tau->psitype = NULL;

    /*  Set the error_occurred member to false and the error_message to NULL. *
     *  If no errors occur during processing, these variables will remain     *
     *  unchanged. Check them throughout to ensure no illegal actions happen. */
    tau->error_occurred = rssringoccs_False;
    tau->error_message = NULL;

    /*  Check if the input dlp is NULL.                                       */
    if (dlp == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rInput dlp is NULL. Returning.\n"
        );
        return tau;
    }

    /*  Check if the input dlp has its error_occurred member set to true.     */
    if (dlp->error_occurred)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rInput dlp is NULL. Returning.\n"
        );
        return tau;
    }

    /*  Set the resolution variable to the user-provided input.               */
    tau->res = res;

    /*  Check that the resolution is a legal value.                           */
    if (tau->res <= 0.0)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rInput res is not positive. Returning.\n"
        );
    }

    /*  Set the default processing keywords.                                  */
    tau->bfac     = rssringoccs_True;
    tau->use_norm = rssringoccs_True;
    tau->use_fft  = rssringoccs_False;
    tau->use_fwd  = rssringoccs_False;

    /*  The default window is the modified Kaiser-Bessel with 2.0 alpha.      */
    tau->wtype = rssringoccs_strdup("kbmd20");

    /*  Check that rssringoccs_strdup did not fail.                           */
    if (tau->wtype == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rrssringoccs_strdup failed to set tau->wtype. Returning.\n"
        );
        return tau;
    }

    /*  Set normeq to the value corresponding to "kbmd20".                    */
    tau->normeq = KBMD20NormEQ;

    /*  Set the window function to kbmd20. This can be changed later.         */
    tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_2_0;

    /*  And the default psitype is Fresnel processing using quartic           *
     *  Legendre polynomials.                                                 */
    tau->psitype = rssringoccs_strdup("fresnel4");

    /*  Check that rssringoccs_strdup did not fail.                           */
    if (tau->psitype == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rrssringoccs_strdup failed to set tau->psitype. Returning.\n"
        );

        return tau;
    }

    /*  Set the order and psinum corresponding to "fresnel4".                 */
    tau->order  = 4;
    tau->psinum = rssringoccs_DR_Legendre;

    /*  The default sigma (Allen Deviation) is the Cassini value.             */
    tau->sigma = 2.0e-13;

    /*  Set the default values for eccentricity and peripase. We always       *
     *  assume circular orbits. The user must specify otherwise.              */
    tau->ecc  = 0.0;
    tau->peri = 0.0;

    /*  By default we do not perturb psi by polynomials. Set these to 0.      */
    tau->perturb[0] = 0.0;
    tau->perturb[1] = 0.0;
    tau->perturb[2] = 0.0;
    tau->perturb[3] = 0.0;
    tau->perturb[4] = 0.0;

    /*  The default range request is "all" which has numerical values of      *
     *  1.0 and 400,000. This is more than enough for Saturn and Uranus       *
     *  occultations. The actual routines will trim this down to the          *
     *  available data and store these in the rng_list member later.          */
    tau->rng_req[0] = 1.0;
    tau->rng_req[1] = 4.0e5;

    /*  By default, we do not use the interpolation methods defined in MTR86. */
    tau->interp = 0;

    /**************************************************************************
     *  Grab the data from the DLP and compute some extra variables. This     *
     *  function computes the following for tau:                              *
     *      dx_km                                                             *
     *      arr_size                                                          *
     *      rho_km_vals                                                       *
     *      phi_rad_vals                                                      *
     *      f_sky_hz_vals                                                     *
     *      rho_dot_kms_vals                                                  *
     *      raw_tau_threshold_vals                                            *
     *      B_rad_vals                                                        *
     *      D_km_vals                                                         *
     *      t_oet_spm_vals                                                    *
     *      t_ret_spm_vals                                                    *
     *      t_set_spm_vals                                                    *
     *      rho_corr_pole_km_vals                                             *
     *      rho_corr_timing_km_vals                                           *
     *      phi_rl_rad_vals                                                   *
     *      p_norm_vals                                                       *
     *      phase_rad_vals                                                    *
     *  This function also checks that the arrays have valid entries.         *
     **************************************************************************/
    rssringoccs_Copy_DLP_Data_To_Tau(dlp, tau);
    return tau;
}
/*  rssringoccs_Create_TAUObj.                                                */
