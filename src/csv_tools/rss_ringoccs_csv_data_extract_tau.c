
/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void
rssringoccs_CSVData_Extract_Tau(rssringoccs_CSVData *csv, const char *tau_file)
{
    /*  Buffer for an error message, should an error occur.                   */
    char err_mes[1024];

    /*  If the input pointer is NULL, there is nothing to do.                 */
    if (!csv)
        return;

    /*  Similarly if an error occurred. Abort the computation.                */
    if (csv->error_occurred)
        return;

    /*  Extract the data from the GEO.TAB file.                               */
    csv->tau = rssringoccs_TauCSV_Extract(tau_file, csv->use_deprecated);

    /*  Check for errors.                                                     */
    if (!csv->tau)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Extract_Tau\n\n"
            "rssringoccs_TauCSV_Extract returned NULL. Aborting.\n"
        );

        return;
    }

    /*  Make sure the extraction was successful.                              */
    if (csv->tau->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous tau message.      */
        if (csv->tau->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_Tau\n\n"
                "rssringoccs_TauCSV_Extract returned with error.\n"
                "rssringoccs_TauCSV_Extract set the following message:\n\n%s",
                csv->tau->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_Tau\n\n"
                "rssringoccs_TauCSV_Extract returned with error.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_CSVData_Destroy_Members(csv);
        return;
    }

    /*  If the CSV is empty there is something wrong with the tau string.     */
    if (csv->tau->n_elements == 0)
    {
        csv->error_occurred = tmpl_True;

        /*  Zero elements should have been treated as an error. Check.        */
        if (csv->tau->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_Tau\n\n"
                "rssringoccs_TauCSV_Extract returned an empty struct.\n"
                "rssringoccs_TauCSV_Extract set the following message:\n\n%s",
                csv->tau->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  Otherwise, treat it as an error now.                              */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_Tau\n\n"
                "rssringoccs_TauCSV_Extract returned an empty struct.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_CSVData_Destroy_Members(csv);
        return;
    }
}
