#ifndef RSS_RINGOCCS_HISTORY_H
#define RSS_RINGOCCS_HISTORY_H

typedef struct rssringoccs_History_Def {
    const char *rss_ringoccs_version;
    const char *libtmpl_version;
    const char *c_version;
    const char *user_name;
    const char *host_name;
    const char *run_date;
    const char *operating_system;

    /*  To avoid having to malloc this array, it is statically sized. The     *
     *  strings come in pairs "parameter_name: the_parameter", and the array  *
     *  ends in a NULL terminator. So you can have (17 - 1) / 2 = 8 different *
     *  arguments, and similarly 8 different keywords. For the rss_ringoccs   *
     *  structs that have history data, this is always sufficient.            */
    const char *input_vars[17];
    const char *input_kwds[17];
} rssringoccs_History;

extern char *
rssringoccs_Date_to_Rev(unsigned int year, unsigned int doy);

extern const char *rssringoccs_Get_Version(void);

extern void
rssringoccs_History_Print(const rssringoccs_History * const history);

extern void
rssringoccs_History_Create_Base(rssringoccs_History *history);

#endif
