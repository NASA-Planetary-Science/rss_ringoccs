/*  NULL pointers are given here.                                             */
#include <stddef.h>

/*  Optical and math functions provided by this library.                      */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for computing several Tau variables from DLP data.               */
void
rssringoccs_Tau_Compute_Data_From_DLP_Members(rssringoccs_TAUObj *tau,
                                              const rssringoccs_DLPObj *dlp)
{

    /*  Constant for zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  Variable for indexing over the data.                                  */
    size_t n;

    /*  Variable for the wavelength of the wave in the DLP object.            */
    double lambda_sky;

    /*  If the tau object is NULL there is nothing to be done.                */
    if (!tau)
        return;

    /*  Similarly if an error occurred before this function was called.       */
    if (tau->error_occurred)
        return;

    /*  The DLP object should not be NULL. Check for this.                    */
    if (!dlp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Data_From_DLP_Members\n\n"
            "\rInput DLP object is NULL. Aborting.\n\n"
        );

        return;
    }

    /*  If the input DLP had an error occur previously, treat this as an      *
     *  error. Store an error message in the Tau object.                      */
    if (dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Data_From_DLP_Members\n\n"
            "\rInput DLP object has error_occurred = True.\n\n"
        );

        return;
    }

    /*  The DLP and Tau object should have the same number of elements        *
     *  allocated for each of their members. Check for this.                  */
    if (dlp->arr_size != tau->arr_size)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Compute_Data_From_DLP_Members\n\n"
            "\rDLP array size is not equal to Tau array size.\n\n"
        );

        return;
    }

    /*  Tau has complex data, DLP has real data. Compute the all the needed   *
     *  variables for the Tau object using this DLP data.                     */
    for (n = zero; n < tau->arr_size; ++n)
    {
        /*  Compute the complex amplitude, T_hat_vals.                        */
        tau->T_in[n] = tmpl_CDouble_Polard(
            tmpl_Double_Sqrt(dlp->p_norm_vals[n]), -dlp->phase_deg_vals[n]
        );

        /*  Compute the wavelength lambda.                                    */
        lambda_sky = tmpl_Double_Frequency_To_Wavelength(dlp->f_sky_hz_vals[n]);

        /*  Use the wagelength to compute the wavenumber.                     */
        tau->k_vals[n] = tmpl_Double_Wavelength_To_Wavenumber(lambda_sky);

        /*  And finally, compute the Fresnel scale.                           */
        tau->F_km_vals[n] = tmpl_Double_Cyl_Fresnel_Scale_Deg(
            lambda_sky,             /*  Wavelength.                           */
            tau->D_km_vals[n],      /*  Spacecraft-to-Ring-Plane Distance.    */
            tau->phi_deg_vals[n],   /*  Ring azimuth angle (degrees).         */
            tau->B_deg_vals[n]      /*  Ring opening angle (degrees).         */
        );
    }
}
/*  End of rssringoccs_Tau_Compute_Data_From_DLP_Members.                     */
