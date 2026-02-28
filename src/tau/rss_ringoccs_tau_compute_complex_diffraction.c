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
 *                rss_ringoccs_tau_compute_complex_diffraction                *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Computes the complex diffracted transmittance for a Tau object.       *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Compute_Complex_Diffraction                           *
 *  Purpose:                                                                  *
 *      Computes the complex diffracted transmittance using the data in a DLP.*
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * TMPL_RESTRICT const):                       *
 *          The Tau object.                                                   *
 *      dlp (const rssringoccs_DLPObj * TMPL_RESTRICT const):                 *
 *          The DLP object.                                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_complex.h:                                                       *
 *          tmpl_CDouble_Polard:                                              *
 *              Creates a complex number from polar form (in degrees).        *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Sqrt:                                                 *
 *              Computes the square root of a real number.                    *
 *  Method:                                                                   *
 *      The DLP object contains the raw power and the phase. The complex      *
 *      diffracted transmittance can be computed from this:                   *
 *                                                                            *
 *          ^     _                                                           *
 *          T = \/p exp(-i theta)                                             *
 *                                                                            *
 *      We compute the square root of the power, and then create T_hat using  *
 *      the polar form of a complex number.                                   *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers before trying to access    *
 *          data. The error_occurred Boolean is set to true if one of the     *
 *          required variables is NULL.                                       *
 *                                                                            *
 *      2.) Both the tau and dlp pointers are declared with TMPL_RESTRICT. On *
 *          compilers supporting the C99 standard, this expands to "restrict" *
 *          meaning tau and dlp must point to different objects. This should  *
 *          be true regardless in order to properly use this function.        *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. T_hat is defined in this paper.           *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans (True, False).                     *
 *  2.) tmpl_complex.h:                                                       *
 *          Header file providing complex numbers and functions.              *
 *  3.) tmpl_math.h:                                                          *
 *          Header file providing the square root function.                   *
 *  4.) rss_ringoccs_tau.h:                                                   *
 *          Header providing the TAU and DLP typedefs, and function prototype.*
 *  5.) stddef.h:                                                             *
 *          Standard library header providing the size_t typedef.             *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       November 28, 2023                                             *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2025/04/11: Ryan Maguire                                                  *
 *      Removed call to tmpl_String_Duplicate, error_message is now a const   *
 *      char pointer, no need to call free or malloc.                         *
 *  2026/02/21: Ryan Maguire                                                  *
 *      Added checks for NULL pointers in the geometry, moved the T_in        *
 *      computation to its own function, changed function name.               *
 *  2026/02/27: Ryan Maguire                                                  *
 *      Changed function signature, the Tau object now has a pointer to a DLP *
 *      object as a member, we do not need DLP as a function parameter.       *
 ******************************************************************************/

/*  Booleans (True / False) provided here.                                    */
#include <libtmpl/include/tmpl_bool.h>

/*  Complex numbers and functions given here.                                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Square root function found here.                                          */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Function for computing the diffracted transmittance from a DLP object.    */
void rssringoccs_Tau_Compute_Complex_Diffraction(rssringoccs_TAUObj * const tau)
{
    /*  Variable for indexing over the data.                                  */
    size_t n;

    /*  If the tau object is NULL there is nothing to be done.                */
    if (!tau)
        return;

    /*  Similarly if an error occurred before this function was called.       */
    if (tau->error_occurred)
        return;

    /*  The DLP object should not be NULL. Check for this.                    */
    if (!tau->dlp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Complex_Diffraction\n\n"
            "\rtau->dlp is NULL.\n\n";

        return;
    }

    /*  If the input DLP had an error occur previously, treat this as an      *
     *  error. Store an error message in the Tau object.                      */
    if (tau->dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Complex_Diffraction\n\n"
            "\rtau->dlp has error_occurred = True.\n\n";

        return;
    }

    /*  The complex diffracted transmittance, T_in, needs the diffracted      *
     *  power and phase from the DLP. Make sure these are not NULL.           */
    if (!tau->dlp->p_norm_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Complex_Diffraction\n\n"
            "\rtau->dlp->p_norm_vals is NULL.\n\n";

        return;
    }

    /*  Same check for the phase angle.                                       */
    if (!tau->dlp->phase_deg_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Complex_Diffraction\n\n"
            "\rtau->dlp->phase_deg_vals is NULL.\n\n";

        return;
    }

    /*  We will be writing to the complex input array (T_in). This is set to  *
     *  NULL when tau is created and then allocated memory later. Make sure   *
     *  it is not NULL before trying to write to it.                          */
    if (!tau->T_in)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Complex_Diffraction\n\n"
            "\rtau->T_in is NULL, cannot write data.\n\n";

        return;
    }

    /*  Loop through and compute the transmittance from the data in dlp.      */
    for (n = 0; n < tau->dlp->arr_size; ++n)
    {
        /*  Extract the power and the phase. By convention, the phase angle   *
         *  for the Tau object (used for diffraction correction) is negative  *
         *  the phase angle stored in the DLP. This ensures that the Fresnel  *
         *  transforms use the correct sign.                                  */
        const double power = tau->dlp->p_norm_vals[n];
        const double phase = -tau->dlp->phase_deg_vals[n];

        /*  power = | T |^2, hence | T | = sqrt(power). Compute.              */
        const double modulus = tmpl_Double_Sqrt(power);

        /*  Compute the complex amplitude, T_hat.                             */
        tau->T_in[n] = tmpl_CDouble_Polard(modulus, phase);
    }
}
/*  End of rssringoccs_Tau_Compute_Complex_Diffraction.                       */
