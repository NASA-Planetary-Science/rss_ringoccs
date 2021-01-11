#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>


/*  Macro for checking if certain pointers in the tau object are NULL.        *
 *  Several of the inputs should NOT be NULL when calling this function. Note *
 *  that since this macro ends with braces, we do not need a semi-colon at    *
 *  the end of a line after calling it.                                       */
#define CHECK_OLD_TAU_MEMBERS(var)                                             \
    if (tau->var == NULL)                                                      \
    {                                                                          \
        tau->error_occurred = rssringoccs_True;                                \
        tau->error_message = rssringoccs_strdup(                               \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Compute_Vars\n\n"                             \
            "\rInput tau has "#var" set to NULL. Returning.\n\n"               \
        );                                                                     \
        return;                                                                \
    }
/*  End of CHECK_OLD_TAU_MEMBERS macro.                                       */

/*  Macro for checking if the members we are computing are set to NULL. At    *
 *  the start of this function the new members SHOULD be NULL since they are  *
 *  set to NULL by rssringoccs_Create_TAUObj and since no other functions     *
 *  should alter them before now.                                             */
#define CHECK_NEW_TAU_MEMBERS(var)                                             \
    if (tau->var != NULL)                                                      \
    {                                                                          \
        tau->error_occurred = rssringoccs_True;                                \
        tau->error_message = rssringoccs_strdup(                               \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Compute_Vars\n\n"                             \
            "\rInput tau->"#var" is not NULL. It is likely you've already\n"   \
            "\rcalled this function and computed tau->T_in. Returning.\n\n"    \
        );                                                                     \
        return;                                                                \
    }
/*  End of CHECK_NEW_TAU_MEMBERS macro.                                       */

/*  Use this macro to save on repetitive code. It attempts to allocate memory *
 *  for a member of a tau object and then checks if malloc failed.            */
#define MALLOC_TAU_MEMBER(var)                                                 \
    /*  Allocate memory for the variable.                                    */\
    tau->var = malloc(sizeof(*tau->var) * tau->arr_size);                      \
                                                                               \
    /*  Check if malloc failed.                                              */\
    if (tau->var == NULL)                                                      \
    {                                                                          \
        tau->error_occurred = rssringoccs_True;                                \
        tau->error_message = rssringoccs_strdup(                               \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"                         \
            "\rMalloc failed and returned NULL for "#var". Returning.\n\n"     \
        );                                                                     \
        return;                                                                \
    }
/*  End of the MALLOC_TAU_MEMBER macro.                                       */

/*  Function for computing T_hat, F, and k for a tau object.                  */
void rssringoccs_Tau_Compute_Vars(rssringoccs_TAUObj *tau)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
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
    CHECK_OLD_TAU_MEMBERS(p_norm_vals)
    CHECK_OLD_TAU_MEMBERS(phase_rad_vals)
    CHECK_OLD_TAU_MEMBERS(f_sky_hz_vals)
    CHECK_OLD_TAU_MEMBERS(D_km_vals)
    CHECK_OLD_TAU_MEMBERS(phi_rad_vals)
    CHECK_OLD_TAU_MEMBERS(B_rad_vals)

    /*  Check that the variables we're computing have not be malloc'd or      *
     *  computed already.                                                     */
    CHECK_NEW_TAU_MEMBERS(T_in)
    CHECK_NEW_TAU_MEMBERS(F_km_vals)
    CHECK_NEW_TAU_MEMBERS(k_vals)

    /*  Allocate memory for the complex transmittance T_in, the Fresnel scale *
     *  F_km_vals, and the wavenumber k_vals.                                 */
    MALLOC_TAU_MEMBER(T_in)
    MALLOC_TAU_MEMBER(F_km_vals)
    MALLOC_TAU_MEMBER(k_vals)

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
