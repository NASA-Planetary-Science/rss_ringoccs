#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_calendar.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_utility.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>

void
rssringoccs_GeoCSV_Write_History(rssringoccs_GeoCSV *geo, const char *filename)
{
    if (!geo)
        return;

    if (geo->error_occurred)
        return;

    geo->history = malloc(sizeof(*geo->history));

    if (!geo->history)
    {
        geo->error_occurred = tmpl_True;
        geo->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_GeoCSV_Write_History\n\n"
            "malloc failed to allocate memory for history.\n";

        return;
    }

    rssringoccs_History_Create_Base(geo->history);

    geo->history->input_vars[0] = "filename";
    geo->history->input_vars[1] = filename;
    geo->history->input_vars[2] = NULL;
    geo->history->input_kwds[0] = "use_deprecated";

    if (geo->use_deprecated)
        geo->history->input_kwds[1] = "True";
    else
        geo->history->input_kwds[1] = "False";

    geo->history->input_kwds[2] = NULL;
}
