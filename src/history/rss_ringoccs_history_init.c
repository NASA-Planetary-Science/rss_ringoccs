#include <libtmpl/include/helper/tmpl_array_size.h>
#include <libtmpl/include/tmpl_calendar.h>
#include <libtmpl/include/tmpl_utility.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <stddef.h>

#ifdef __STDC_VERSION__
#define RSSRINGOCCS_TO_STRING(x) #x
#define RSSRINGOCCS_MAKE_STRING(x) RSSRINGOCCS_TO_STRING(x)
#endif

void
rssringoccs_History_Init(rssringoccs_History * const history)
{
    size_t n;

    if (!history)
        return;

    history->rss_ringoccs_version = rssringoccs_Version();
    history->libtmpl_version = tmpl_Version();

#ifndef __STDC_VERSION__
    history->c_version = "Unknown";
#else
    history->c_version = RSSRINGOCCS_MAKE_STRING(__STDC_VERSION__);
#endif

    history->user_name = tmpl_User_Name();
    history->host_name = tmpl_Host_Name();
    history->run_date = tmpl_Local_Calendar_Date_And_Time();
    history->operating_system = tmpl_Operating_System();

    for (n = 0; n < TMPL_ARRAY_SIZE(history->input_vars); ++n)
        history->input_vars[n] = NULL;

    for (n = 0; n < TMPL_ARRAY_SIZE(history->input_kwds); ++n)
        history->input_kwds[n] = NULL;
}

#undef RSSRINGOCCS_TO_STRING
#undef RSSRINGOCCS_MAKE_STRING
