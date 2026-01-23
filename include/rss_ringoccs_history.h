#ifndef RSS_RINGOCCS_HISTORY_H
#define RSS_RINGOCCS_HISTORY_H

#include <rss_ringoccs/include/types/rss_ringoccs_history.h>

extern const char *
rssringoccs_Date_to_Rev(const unsigned int year, const unsigned int doy);

extern const char *rssringoccs_Get_Version(void);

extern void
rssringoccs_History_Print(const rssringoccs_History * const history);

extern void
rssringoccs_History_Create_Base(rssringoccs_History *history);

#endif
