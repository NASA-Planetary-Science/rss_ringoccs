
/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void
rssringoccs_CSVData_Extract_Geo(rssringoccs_CSVData *csv, const char *geo_file)
{
    /*  If the input pointer is NULL, there is nothing to do.                 */
    if (!csv)
        return;

    /*  Similarly if an error occurred. Abort the computation.                */
    if (csv->error_occurred)
        return;

    /*  Extract the data from the GEO.TAB file.                               */
    csv->geo = rssringoccs_GeoCSV_Extract(geo_file, csv->use_deprecated);

    /*  Check for errors.                                                     */
    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Extract_Geo\n\n"
            "rssringoccs_GeoCSV_Extract returned NULL.\n\n";

        return;
    }

    /*  Make sure the extraction was successful.                              */
    if (csv->geo->error_occurred)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Extract_Geo\n\n"
            "rssringoccs_GeoCSV_Extract returned with error.\n\n";

        /*  Free all data and abort.                                          */
        rssringoccs_CSVData_Destroy_Members(csv);
        return;
    }

    /*  If the CSV is empty there is something wrong with the geo string.     */
    if (csv->geo->n_elements == 0)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Extract_Geo\n\n"
                "rssringoccs_GeoCSV_Extract returned an empty struct.\n\n";

        /*  Free all data and abort.                                          */
        rssringoccs_CSVData_Destroy_Members(csv);
        return;
    }
}
