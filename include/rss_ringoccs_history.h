

#ifndef RSS_RINGOCCS_HISTORY_H
#define RSS_RINGOCCS_HISTORY_H

typedef struct rssringoccs_HistoryObj_Def {
    char *rssocc_version;
    char *user_name;
    char *host_name;
    char *run_date;
    char *os;
    char **inpute_vars;
    char **input_kwds;
} rssringoccs_HistoryObj;

extern char *
rssringoccs_Date_to_Rev(unsigned int year, unsigned int doy);

#endif
