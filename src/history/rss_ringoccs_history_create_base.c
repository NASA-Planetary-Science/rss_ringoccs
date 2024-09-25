#include <libtmpl/include/tmpl_calendar.h>
#include <libtmpl/include/tmpl_utility.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <stdlib.h>

#define RSSRINGOCCS_TO_STRING(x) #x
#define RSSRINGOCCS_MAKE_STRING(x) RSSRINGOCCS_TO_STRING(x)

void
rssringoccs_History_Create_Base(rssringoccs_History *history)
{
    if (!history)
        return;

    history->rss_ringoccs_version = rssringoccs_Get_Version();
#ifndef __STDC_VERSION__
    history->c_version = "Unknown";
#else
    history->c_version = RSSRINGOCCS_MAKE_STRING(__STDC_VERSION__);
#endif
    history->libtmpl_version = tmpl_Version();
    history->host_name = tmpl_Host_Name();
    history->user_name = tmpl_User_Name();
    history->run_date = tmpl_Local_Calendar_Date_And_Time();
    history->operating_system = tmpl_Operating_System();
}

#undef RSSRINGOCCS_TO_STRING
#undef RSSRINGOCCS_MAKE_STRING
