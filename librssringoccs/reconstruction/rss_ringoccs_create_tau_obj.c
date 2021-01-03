#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

/*  Function for allocating memory for a Tau object and setting the default   *
 *  values for all of the keywords.                                           */
rssringoccs_TAUObj *rssringoccs_Create_TAUObj(void)
{
    rssringoccs_TAUObj *tau;
    tau = malloc(sizeof(*tau));

    if (tau == NULL)
        return tau;

    /*  Initialize all pointers within the tau struct to NULL.                */
    tau->T_in = NULL;
    tau->T_out = NULL;
    tau->rho_km_vals = NULL;
    tau->F_km_vals = NULL;
    tau->phi_rad_vals = NULL;
    tau->kd_vals = NULL;
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

    /*  Set the default processing keywords.                                  */
    tau->bfac     = rssringoccs_True;
    tau->use_norm = rssringoccs_True;
    tau->use_fft  = rssringoccs_False;
    tau->use_fwd  = rssringoccs_False;

    /*  If processing proceeds without error, this will remain at false.      */
    tau->error_occurred = rssringoccs_False;
    tau->error_message = NULL;

    /*  The default window is the modified Kaiser-Bessel with 2.0 alpha.      */
    tau->wtype = rssringoccs_strdup("kbmd20");

    /*  Check that rssringoccs_strdup did not fail.                           */
    if (tau->wtype == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Create_TAUObj\n\n"
            "\rrssringoccs_strdup failed to set tau->wtype. Returning.\n"
        );

        /*  Set psitype to NULL before returning so we don't try to free a    *
         *  pointer that wasn't malloc'd.                                     */
        tau->psitype = NULL;
        return tau;
    }

    /*  And the default psitype is Fresnel processing using quartic           *
     *  Legendre polynomials.                                                 */
    tau->psitype = rssringoccs_strdup("fresnel4");
    tau->order   = 4;
    tau->psinum  = rssringoccs_DR_Legendre;

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

    /*  Set the default values for eccentricity and peripase. We always       *
     *  assume circular orbits. The user must specify otherwise.              */
    tau->ecc  = 0.0;
    tau->peri = 0.0;

    /*  By default we do not perturb psi by polynomials. Set coefficient to 0.*/
    tau->perturb[0] = 0.0;
    tau->perturb[1] = 0.0;
    tau->perturb[2] = 0.0;
    tau->perturb[3] = 0.0;
    tau->perturb[4] = 0.0;

    /*  Initialize dx to -1. If processing proceeds and dx is still -1, an    *
     *  error occured and the reconstruction will be aborted.                 */
    tau->dx = -1.0;

    /*  By default, we do not use the interpolation methods defined in MTR86. */
    tau->interp = 0;

    return tau;
}
/*  rssringoccs_Create_TAUObj.                                                */
