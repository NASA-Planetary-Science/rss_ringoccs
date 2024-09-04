#include <stdlib.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Macro function for safely allocating memory for the variables. This       *
 *  checks if malloc fails, and does not simply assume it passed.             */
#define MALLOC_CSV_VAR(var)                                                    \
    csv->var = malloc(sizeof(*csv->var)*csv->n_elements);                      \
    if (csv->var == NULL)                                                      \
    {                                                                          \
        csv->error_occurred = tmpl_True;                                       \
        csv->error_message = tmpl_String_Duplicate(                            \
            "Error Encountered: rss_ringoccs\n"                                \
            "\trssringoccs_CSVData_Malloc\n\n"                                 \
            "Malloc returned NULL for csv member. Aborting.\n"                 \
        );                                                                     \
                                                                               \
        rssringoccs_CSVData_Destroy_Members(csv);                              \
        return;                                                                \
    }

void rssringoccs_CSVData_Malloc(rssringoccs_CSVData *csv)
{
    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (csv->n_elements == (size_t)0)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Malloc\n\n"
            "n_elements is zero, nothing to malloc. Aborting.\n"
        );

        return;
    }

    /*  The MALLOC_CSV_VAR macro contains an if-then statement with           *
     *  braces {} hence there is no need for a semi-colon at the end.         */
    MALLOC_CSV_VAR(D_km_vals)
    MALLOC_CSV_VAR(f_sky_hz_vals)
    MALLOC_CSV_VAR(rho_dot_kms_vals)
    MALLOC_CSV_VAR(rx_km_vals)
    MALLOC_CSV_VAR(ry_km_vals)
    MALLOC_CSV_VAR(rz_km_vals)

    /*  If Tau data is to be extracted, reserve memory for the variables.     */
    if (csv->tau)
    {
        MALLOC_CSV_VAR(tau_phase)
        MALLOC_CSV_VAR(tau_power)
        MALLOC_CSV_VAR(tau_vals)
    }
}
