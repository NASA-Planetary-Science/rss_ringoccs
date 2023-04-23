

#ifndef RSS_RINGOCCS_HISTORY_H
#define RSS_RINGOCCS_HISTORY_H

#include "librssringoccs_exports.h"

typedef struct rssringoccs_HistoryObj_Def {
    char *rssocc_version;
    char *user_name;
    char *host_name;
    char *run_date;
    char *os;
    char **inpute_vars;
    char **input_kwds;
} rssringoccs_HistoryObj;

RSS_RINGOCCS_EXPORT extern char *
rssringoccs_Date_to_Rev(unsigned int year, unsigned int doy);

#endif
