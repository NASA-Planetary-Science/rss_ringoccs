
/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void
rssringoccs_UranusCSVData_Extract_Geo(rssringoccs_UranusCSVData *csv,
                                      const char *geo_file)
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
    csv->geo = rssringoccs_GeoCSV_Extract(geo_file, tmpl_False);

    /*  Check for errors.                                                     */
    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Extract_Geo\n\n"
            "rssringoccs_GeoCSV_Extract returned NULL. Aborting.\n"
        );

        return;
    }

    /*  Make sure the extraction was successful.                              */
    if (csv->geo->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous geo message.      */
        if (csv->geo->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_Geo\n\n"
                "rssringoccs_GeoCSV_Extract returned with error.\n"
                "rssringoccs_GeoCSV_Extract set the following message:\n\n%s",
                csv->geo->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_Geo\n\n"
                "rssringoccs_GeoCSV_Extract returned with error.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_UranusCSVData_Destroy_Members(csv);
        return;
    }

    /*  If the CSV is empty there is something wrong with the geo string.     */
    if (csv->geo->n_elements == 0)
    {
        csv->error_occurred = tmpl_True;

        /*  Zero elements should have been treated as an error. Check.        */
        if (csv->geo->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_Geo\n\n"
                "rssringoccs_GeoCSV_Extract returned an empty struct.\n"
                "rssringoccs_GeoCSV_Extract set the following message:\n\n%s",
                csv->geo->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  Otherwise, treat it as an error now.                              */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Extract_Geo\n\n"
                "rssringoccs_GeoCSV_Extract returned an empty struct.\n"
            );
        }

        /*  Free all data and abort.                                          */
        rssringoccs_UranusCSVData_Destroy_Members(csv);
        return;
    }
}
