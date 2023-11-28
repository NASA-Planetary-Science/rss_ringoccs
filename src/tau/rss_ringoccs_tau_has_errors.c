/*  NULL pointers are given here.                                             */
#include <stddef.h>

/*  Optical and math functions provided by this library.                      */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  As a side note, the expression #var is used in the following macros. This *
 *  takes the input of the macro, which we're calling var for all three, and  *
 *  treats it as a string literal.                                            */

/*  Use this macro to save on repetitive code. It is for checking that all of *
 *  of the values of a given member in the tau object are non-negative.       */
#define TAU_CHECK_NON_NEGATIVE(var)                                            \
                                                                               \
    /*  Use tmpl_Double_Array_Min to compute the minimum value and check if  */\
    /*  it is negative. Return error if it is.                               */\
    if (tmpl_Double_Array_Min(tau->var, tau->arr_size) < 0.0)                  \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_String_Duplicate(                            \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Has_Errors\n\n"                               \
            "\r"#var" has negative valued entries. Returning.\n"               \
        );                                                                     \
                                                                               \
        return tmpl_True;                                                      \
    }
/*  End of the TAU_CHECK_NON_NEGATIVE macro.                                  */

/*  Use this macro to save on repetitive code. It is for checking that all of *
 *  of the values of a given member in the tau object fall within [-2pi, 2pi].*
 *  Note, it implicitly has min and max defined. These are declared at the    *
 *  top of the rssringoccs_Copy_DLP_Data_To_Tau function.                     */
#define TAU_CHECK_360(var)                                                     \
                                                                               \
    /*  Compute the minimum and maximum of var.                              */\
    tmpl_Double_Array_MinMax(tau->var, tau->arr_size, &min, &max);             \
                                                                               \
    /*  Check if var falls within the interval [-2pi, 2pi].                  */\
    if ((min < -360.0) || (max > 360.0))                                       \
    {                                                                          \
        tau->error_occurred = tmpl_True;                                       \
        tau->error_message = tmpl_String_Duplicate(                            \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\trssringoccs_Tau_Has_Errors\n\n"                               \
            "\r"#var" has values outside of [-360, 360]. Returning.\n"         \
        );                                                                     \
                                                                               \
        return tmpl_True;                                                      \
    }
/*  End of the TAU_CHECK_TWO_PI macro.                                        */

/*  Function for checking Tau parameters for possible errors.                 */
tmpl_Bool rssringoccs_Tau_Has_Errors(rssringoccs_TAUObj *tau)
{
    /*  Variables for keeping track of the min and max of arrays.             */
    double min, max;

    if (!tau)
        return tmpl_True;

    if (tau->error_occurred)
        return tmpl_True;

    /*  Check if dx_km is a legal value.                                      */
    if (tau->dx_km == 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Has_Errors\n\n"
            "\rdx_km = rho_km_vals[1] - rho_km_vals[0] = 0. Returning.\n"
        );

        return tmpl_True;
    }

    /*  dx_km may be negative if this is an ingress occultation. To check if  *
     *  res is a legal value, compare it with twice the absolute value of     *
     *  dx_km. To avoid floating round-off error (which has happened to the   *
     *  Cassini team, hence this edit) set the value to 1.99 instead of 2.0.  */
    else if (tau->res < 1.99 * tmpl_Double_Abs(tau->dx_km))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Has_Errors\n\n"
            "\rResolution is less than twice the sample space.\n"
            "\rThis will result in an inaccurate reconstruction. Returning.\n"
        );

        return tmpl_True;
    }

    /*  The following members of tau should be non-negative for all entries.  *
     *  The __TAU_CHECK_NON_NEGATIVE__ contains an if statement and ends with *
     *  braces {}. Hence we do not need a semi-colon at the end.              */
    TAU_CHECK_NON_NEGATIVE(rho_km_vals)
    TAU_CHECK_NON_NEGATIVE(D_km_vals)

    /*  Check that the following variables for angles fall within [-2pi, 2pi].*
     *  Like the other two macros, the __TAU_CHECK_TWO_PI__ macro ends with   *
     *  braces {} so we do not need a semi-colon at the end of these lines.   */
    TAU_CHECK_360(B_deg_vals)
    TAU_CHECK_360(phi_deg_vals)

    /*  If we made it this far, the Tau object should be good to go.          */
    return tmpl_False;
}
/*  End of rssringoccs_Tau_Has_Errors.                                        */

#undef TAU_CHECK_NON_NEGATIVE
#undef TAU_CHECK_TWO_PI
