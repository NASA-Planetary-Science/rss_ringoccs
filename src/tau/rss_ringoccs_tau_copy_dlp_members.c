/*  NULL pointers are given here.                                             */
#include <stddef.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

void
rssringoccs_Tau_Copy_DLP_Members(rssringoccs_TAUObj *tau,
                                 const rssringoccs_DLPObj *dlp)
{
    /*  Constant for zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  Variable for indexing over the data.                                  */
    size_t n;

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
            "\r\trssringoccs_Tau_Copy_DLP_Members\n\n"
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
            "\r\trssringoccs_Tau_Copy_DLP_Members\n\n"
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
            "\r\trssringoccs_Tau_Copy_DLP_Members\n\n"
            "\rDLP array size is not equal to Tau array size.\n\n"
        );

        return;
    }

    /*  Copy the data from the DLP object to the Tau object.                  */
    for (n = zero; n < tau->arr_size; ++n)
    {
        tau->rho_km_vals[n] = dlp->rho_km_vals[n];
        tau->phi_deg_vals[n] = dlp->phi_deg_vals[n];
        tau->B_deg_vals[n] = dlp->B_deg_vals[n];
        tau->D_km_vals[n] = dlp->D_km_vals[n];
        tau->rho_dot_kms_vals[n] = dlp->rho_dot_kms_vals[n];
        tau->t_oet_spm_vals[n] = dlp->t_oet_spm_vals[n];
        tau->t_ret_spm_vals[n] = dlp->t_ret_spm_vals[n];
        tau->t_set_spm_vals[n] = dlp->t_set_spm_vals[n];
        tau->rho_corr_pole_km_vals[n] = dlp->rho_corr_pole_km_vals[n];
        tau->rho_corr_timing_km_vals[n] = dlp->rho_corr_timing_km_vals[n];
        tau->phi_rl_deg_vals[n] = dlp->phi_rl_deg_vals[n];
        tau->rx_km_vals[n] = dlp->rx_km_vals[n];
        tau->ry_km_vals[n] = dlp->ry_km_vals[n];
        tau->rz_km_vals[n] = dlp->rz_km_vals[n];
    }
}
/*  End of rssringoccs_Tau_Copy_DLP_Members.                                  */
