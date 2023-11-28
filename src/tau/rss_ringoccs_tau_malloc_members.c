/*  NULL pointer and malloc are given here.                                   */
#include <stdlib.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Use this macro to save on repetitive code. It checks if tau->var is NULL, *
 *  attempts to malloc memory for tau->var if it is, and then checks to see   *
 *  if malloc failed.                                                         */
#define MALLOC_TAU_VAR(var)                                                    \
                                                                               \
    /*  Check if the variable is not NULL. It should be at the start.        */\
    if (tau->var != NULL)                                                      \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_String_Duplicate(                            \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Malloc_Members\n\n"                           \
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
        tau->error_message = tmpl_String_Duplicate(                            \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Malloc_Members\n\n"                           \
            "\rMalloc failed and returned NULL for "#var". Returning.\n\n"     \
        );                                                                     \
        return;                                                                \
    }
/*  End of the MALLOC_TAU_VAR macro.                                          */

/*  Function for allocating memory for all of the tau variables.              */
void rssringoccs_Tau_Malloc_Members(rssringoccs_TAUObj *tau)
{
    const size_t zero = (size_t)0;

    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    if (tau->arr_size == zero)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Malloc_Members\n\n"
            "\rInput tau has arr_size = zero. Nothing to allocate.\n\n"
        );

        return;
    }

    /*  The MALLOC_TAU_VAR macro ends with an if statement and so has         *
     *  braces {}. Because of this, we do not need a semi-colon at the end.   *
     *  This macro allocates memory for the members of the tau object and     *
     *  checks for errors.                                                    */
    MALLOC_TAU_VAR(rho_km_vals)
    MALLOC_TAU_VAR(phi_deg_vals)
    MALLOC_TAU_VAR(k_vals)
    MALLOC_TAU_VAR(rho_dot_kms_vals)
    MALLOC_TAU_VAR(B_deg_vals)
    MALLOC_TAU_VAR(D_km_vals)
    MALLOC_TAU_VAR(t_oet_spm_vals)
    MALLOC_TAU_VAR(t_ret_spm_vals)
    MALLOC_TAU_VAR(t_set_spm_vals)
    MALLOC_TAU_VAR(rho_corr_pole_km_vals)
    MALLOC_TAU_VAR(rho_corr_timing_km_vals)
    MALLOC_TAU_VAR(phi_rl_deg_vals)
    MALLOC_TAU_VAR(rx_km_vals)
    MALLOC_TAU_VAR(ry_km_vals)
    MALLOC_TAU_VAR(rz_km_vals)
    MALLOC_TAU_VAR(T_in)
    MALLOC_TAU_VAR(F_km_vals)
}
/*  End of rssringoccs_Tau_Malloc_Members.                                    */

#undef MALLOC_TAU_VAR
