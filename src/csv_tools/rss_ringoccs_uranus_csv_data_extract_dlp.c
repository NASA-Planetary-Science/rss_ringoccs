
/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void
rssringoccs_UranusCSVData_Extract_DLP(rssringoccs_UranusCSVData *csv,
                                      const char *dlp_file)
{
    /*  Buffer for an error message, should an error occur.                   */
    char err_mes[1024];

    /*  If the input pointer is NULL, there is nothing to do.                 */
    if (!csv)
        return;

    /*  Similarly if an error occurred. Abort the computation.                */
    if (csv->error_occurred)
        return;

    /*  Extract the data from the DLP.TAB file.                               */
    csv->dlp = rssringoccs_UranusDLPCSV_Extract(dlp_file, csv->dlp_in_radians);

    /*  Check for errors.                                                     */
    if (!csv->dlp)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Extract_DLP\n\n"
            "rssringoccs_DLPCSV_Extract returned NULL. Aborting.\n"
        );

        return;
    }

    /*  Make sure the extraction was successful.                              */
    if (csv->dlp->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous dlp message.      */
        if (csv->dlp->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_DLP\n\n"
                "rssringoccs_DLPCSV_Extract returned with error.\n"
                "rssringoccs_DLPCSV_Extract set the following message:\n\n%s",
                csv->dlp->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_DLP\n\n"
                "rssringoccs_DLPCSV_Extract returned with error.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_UranusCSVData_Destroy_Members(csv);
        return;
    }

    /*  If the CSV is empty there is something wrong with the dlp string.     */
    if (csv->dlp->n_elements == 0)
    {
        csv->error_occurred = tmpl_True;

        /*  Zero elements should have been treated as an error. Check.        */
        if (csv->dlp->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_DLP\n\n"
                "rssringoccs_DLPCSV_Extract returned an empty struct.\n"
                "rssringoccs_DLPCSV_Extract set the following message:\n\n%s",
                csv->dlp->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  Otherwise, treat it as an error now.                              */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_DLP\n\n"
                "rssringoccs_DLPCSV_Extract returned an empty struct.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_UranusCSVData_Destroy_Members(csv);
        return;
    }

    /*  Grab the number of elements from the DLP CSV. This will be the number *
     *  of elements in the output.                                            */
    csv->n_elements = csv->dlp->n_elements;
}
