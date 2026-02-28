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
 *                   rss_ringoccs_tau_compute_frensel_scale                   *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Computes the Fresnel scale for a Tau object.                          *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Tau_Compute_Fresnel_Scale                                 *
 *  Purpose:                                                                  *
 *      Computes the Fresnel scale for a Tau object using the data in a DLP.  *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The Tau object.                                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_optics.h:                                                        *
 *          tmpl_Double_Frequency_To_Wavelength:                              *
 *              Converts frequency (Hz) to wavelength (km).                   *
 *          tmpl_Double_Wavelength_To_Wavenumber:                             *
 *              Converts wavelength to wavenumber, k = 2pi / lambda.          *
 *      tmpl_cyl_fresnel_optics.h:                                            *
 *          tmpl_Double_Cyl_Fresnel_Scale_Deg:                                *
 *              Computes the Cylindrical Fresnel scale with angles in degrees.*
 *      stdio.h:                                                              *
 *          puts:                                                             *
 *              Prints a string to the screen.                                *
 *  Method:                                                                   *
 *      The DLP object has the (sky) frequency of the wave and the geometry   *
 *      of the occultation observation. From this the wavenumber, k, and the  *
 *      Fresnel scale, F, can be computed. libtmpl provides the tools to do   *
 *      this, we convert frequency to wavelength, and then pass the geometry  *
 *      values and the wavelength to libtmpl.                                 *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers before trying to access    *
 *          data. The error_occurred Boolean is set to true if one of the     *
 *          required variables is NULL.                                       *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. The Fresnel scale is defined here.        *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans (True, False).                     *
 *  2.) tmpl_cyl_fresnel_optics.h:                                            *
 *          Header file providing tools for Fresnel optics.                   *
 *  3.) tmpl_optics.h:                                                        *
 *          Header file providing basic optics routines.                      *
 *  4.) rss_ringoccs_tau.h:                                                   *
 *          Header providing the TAU and DLP typedefs, and function prototype.*
 *  5.) stddef.h:                                                             *
 *          Standard library header providing the size_t typedef.             *
 *  6.) stdio.h:                                                              *
 *          Standard library header file providing the puts function.         *
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

/*  Cylindrical Fresnel optics routines found here.                           */
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>

/*  Wavelength and Wavenumber conversion functions found here.                */
#include <libtmpl/include/tmpl_optics.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Function for computing the Fresnel scale from a given DLP object.         */
void rssringoccs_Tau_Compute_Fresnel_Scale(rssringoccs_TAUObj * const tau)
{
    /*  Variable for indexing over the data.                                  */
    size_t n;

    /*  If the tau object is NULL there is nothing to be done.                */
    if (!tau)
        return;

    /*  Similarly if an error occurred before this function was called.       */
    if (tau->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (tau->dlp->verbose)
        puts("\r\tTAU: Computing Fresnel scale...");

    /*  The DLP object should not be NULL. Check for this.                    */
    if (!tau->dlp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rInput DLP object is NULL.\n\n";

        return;
    }

    /*  If the input DLP had an error occur previously, treat this as an      *
     *  error. Store an error message in the Tau object.                      */
    if (tau->dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rInput DLP object has error_occurred = True.\n\n";

        return;
    }

    /*  The Fresnel scale needs the frequency variable from the DLP. Ensure   *
     *  that it is not NULL before accessing.                                 */
    if (!tau->dlp->f_sky_hz_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rdlp->f_sky_hz_vals is NULL.\n\n";

        return;
    }

    /*  Similarly, the spacecraft-to-ring-intercept-point distance is needed. */
    if (!tau->dlp->D_km_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rdlp->D_km_vals is NULL.\n\n";

        return;
    }

    /*  Next check the azimuth angle, phi.                                    */
    if (!tau->dlp->phi_deg_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rdlp->phi_deg_vals is NULL.\n\n";

        return;
    }

    /*  Lastly, the ring opening angle, B.                                    */
    if (!tau->dlp->B_deg_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rdlp->B_deg_vals is NULL.\n\n";

        return;
    }

    /*  We will be writing to the wavenumber array (k_vals) and the Fresnel   *
     *  scale array (F_km_vals). These are initialized to NULL when tau is    *
     *  created and then malloc'd space later. Make sure they are not NULL.   */
    if (!tau->k_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rtau->k_vals is NULL, cannot write data.\n\n";

        return;
    }

    /*  Same check for the Fresnel scale.                                     */
    if (!tau->F_km_vals)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Fresnel_Scale\n\n"
            "\rtau->F_km_vals is NULL, cannot write data.\n\n";

        return;
    }

    /*  Loop through and compute the Fresnel scale from the data in dlp.      */
    for (n = 0; n < tau->dlp->arr_size; ++n)
    {
        /*  Compute the wavelength lambda.                                    */
        const double frequency = tau->dlp->f_sky_hz_vals[n];
        const double lambda = tmpl_Double_Frequency_To_Wavelength(frequency);

        /*  Use the wavelength to compute the wavenumber.                     */
        tau->k_vals[n] = tmpl_Double_Wavelength_To_Wavenumber(lambda);

        /*  And finally, compute the Fresnel scale.                           */
        tau->F_km_vals[n] = tmpl_Double_Cyl_Fresnel_Scale_Deg(
            lambda,                     /*  Wavelength.                       */
            tau->dlp->D_km_vals[n],     /*  Spacecraft-to-Ring-Plane Distance.*/
            tau->dlp->phi_deg_vals[n],  /*  Ring azimuth angle (degrees).     */
            tau->dlp->B_deg_vals[n]     /*  Ring opening angle (degrees).     */
        );
    }
}
/*  End of rssringoccs_Tau_Compute_Fresnel_Scale.                             */
