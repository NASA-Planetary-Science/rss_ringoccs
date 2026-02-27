/*  NULL pointer and malloc are given here.                                   */
#include <stdlib.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Use this macro to save on repetitive code. It checks if dlp->var is NULL, *
 *  attempts to malloc memory for dlp->var if it is, and then checks to see   *
 *  if malloc failed.                                                         */
#define MALLOC_DLP_VAR(var)                                                    \
    do {                                                                       \
        /*  Check if the variable is not NULL. It should be at the start.    */\
        if (dlp->var)                                                          \
        {                                                                      \
            dlp->error_occurred = tmpl_True;                                   \
            dlp->error_message =                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\trssringoccs_DLP_Malloc_Members\n\n"                       \
                "\r"#var" is not NULL. It is likely you've already set the\n"  \
                "\rdata for this DLP object.\n\n";                             \
            return;                                                            \
        }                                                                      \
                                                                               \
        /*  Allocate memory for the variable.                                */\
        dlp->var = malloc(sizeof(*dlp->var) * dlp->arr_size);                  \
                                                                               \
        /*  Check if malloc failed.                                          */\
        if (!dlp->var)                                                         \
        {                                                                      \
            dlp->error_occurred = tmpl_True;                                   \
            dlp->error_message =                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\trssringoccs_DLP_Malloc_Members\n\n"                       \
                "\rMalloc failed and returned NULL for "#var".\n\n";           \
            return;                                                            \
        }                                                                      \
    } while (0)
/*  End of the MALLOC_DLP_VAR macro.                                          */

/*  Function for allocating memory for all of the dlp variables.              */
void rssringoccs_DLP_Malloc_Members(rssringoccs_DLPObj *dlp)
{
    if (!dlp)
        return;

    if (dlp->error_occurred)
        return;

    if (dlp->arr_size == 0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Malloc_Members\n\n"
            "\rInput dlp has arr_size = 0, nothing to allocate.\n\n";

        return;
    }

    /*  The MALLOC_DLP_VAR macro ends with an if statement and so has         *
     *  braces {}. Because of this, we do not need a semi-colon at the end.   *
     *  This macro allocates memory for the members of the DLP object and     *
     *  checks for errors.                                                    */
    MALLOC_DLP_VAR(rho_km_vals);
    MALLOC_DLP_VAR(phi_deg_vals);
    MALLOC_DLP_VAR(B_deg_vals);
    MALLOC_DLP_VAR(D_km_vals);
    MALLOC_DLP_VAR(f_sky_hz_vals);
    MALLOC_DLP_VAR(rho_dot_kms_vals);
    MALLOC_DLP_VAR(t_oet_spm_vals);
    MALLOC_DLP_VAR(t_ret_spm_vals);
    MALLOC_DLP_VAR(t_set_spm_vals);
    MALLOC_DLP_VAR(rho_corr_pole_km_vals);
    MALLOC_DLP_VAR(rho_corr_timing_km_vals);
    MALLOC_DLP_VAR(phi_rl_deg_vals);
    MALLOC_DLP_VAR(p_norm_vals);
    MALLOC_DLP_VAR(phase_deg_vals);
    MALLOC_DLP_VAR(raw_tau_threshold_vals);
    MALLOC_DLP_VAR(rx_km_vals);
    MALLOC_DLP_VAR(ry_km_vals);
    MALLOC_DLP_VAR(rz_km_vals);
}
/*  End of rssringoccs_DLP_Malloc_Members.                                    */

#undef MALLOC_DLP_VAR
