#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

void rssringoccs_Tau_Compute_Vars(rssringoccs_TAUObj *tau)
{
    unsigned long n;
    double lambda_sky;

    /*  Check if the tau pointer is NULL, returning if it is.                 */
    if (tau == NULL)
        return;

    /*  Check if the tau->error_occurred member has been set to true. If it   *
     *  is, do not perform any computation and return with error.             */
    if (tau->error_occurred)
        return;

    /*  Check if the required pointers in tau have been set. If not, abort.   */
    if (tau->p_norm_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has p_norm_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->phase_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has phase_rad_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->f_sky_hz_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has f_sky_hz_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->D_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has D_km_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->phi_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has phi_rad_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    else if (tau->B_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau has B_rad_vals set to NULL. Returning.\n\n"
        );
        return;
    }
    /*  End of error checking for the pointers in tau.                        */

    /*  Check that the variables we're computing have not be malloc'd or      *
     *  computed already.                                                     */
    if (tau->T_in != NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau->T_in is not NULL. It is likely you've already\n"
            "\rcalled this function and computed tau->T_in. Returning.\n\n"
        );
        return;
    }
    else if (tau->F_km_vals != NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau->F_km_vals is not NULL. It is likely you've already\n"
            "\rcalled this function and computed tau->F_km_vals. Returning.\n\n"
        );
        return;
    }
    else if (tau->k_vals != NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rInput tau->k_vals is not NULL. It is likely you've already\n"
            "\rcalled this function and computed tau->k_vals. Returning.\n\n"
        );
        return;
    }
    /*  End of error check for the variables we're computing.                 */

    /*  Allocate memory for the complex transmittance.                        */
    tau->T_in = malloc(sizeof(*tau->T_in)      * tau->arr_size);

    /*  Check that malloc didn't fail. Return error if it did.                */
    if (tau->T_in == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rMalloc failed and returned NULL for tau->T_in. Returning.\n\n"
        );
        return;
    }

    /*  Allocate memory for the Fresnel scale.                                */
    tau->F_km_vals = malloc(sizeof(*tau->F_km_vals) * tau->arr_size);

    /*  Check that malloc didn't fail. Return error if it did.                */
    if (tau->F_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rMalloc failed and returned NULL for tau->F_km_vals.\n"
            "\rReturning.\n\n"
        );
        return;
    }

    /*  Finally, the wavenumber.                                              */
    tau->k_vals = malloc(sizeof(*tau->k_vals) * tau->arr_size);

    /*  Check that malloc didn't fail. Return error if it did.                */
    if (tau->k_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Vars\n\n"
            "\rMalloc failed and returned NULL for tau->k_vals.\n"
            "\rReturning.\n\n"
        );
        return;
    }

    /*  Loop over the entries of each pointer and compute their values.       */
    for (n = 0; n < tau->arr_size; ++n)
    {
        /*  Compute the complex amplitude, T_hat_vals.                        */
        tau->T_in[n] = rssringoccs_CDouble_Polar(
            rssringoccs_Double_Sqrt(tau->p_norm_vals[n]),
            tau->phase_rad_vals[n]
        );

        /*  Compute the wavelength lambda.                                    */
        lambda_sky =
            rssringoccs_Double_Frequency_To_Wavelength(tau->f_sky_hz_vals[n]);

        /*  Use the wagelength to compute the wavenumber.                     */
        tau->k_vals[n]
            = rssringoccs_Double_Wavelength_To_Wavenumber(lambda_sky);

        /*  And finally, compute the Fresnel scale.                           */
        tau->F_km_vals[n]
            = rssringoccs_Double_Fresnel_Scale(lambda_sky,
                                               tau->D_km_vals[n],
                                               tau->phi_rad_vals[n],
                                               tau->B_rad_vals[n]);
    }
    /*  End of for loop for k, T, and F.                                      */
}
/*  End of rssringoccs_Tau_Compute_Vars.                                      */
