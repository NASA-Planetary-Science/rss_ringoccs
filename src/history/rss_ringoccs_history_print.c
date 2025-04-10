#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <libtmpl/include/helper/tmpl_array_size.h>
#include <stdio.h>

void
rssringoccs_History_Print(const rssringoccs_History * const history)
{
    if (!history)
    {
        puts("None");
        return;
    }

    printf(
        "rss_ringoccs Version: %s\n"
        "libtmpl Version:      %s\n"
        "C Version:            %s\n"
        "Host Name:            %s\n"
        "User Name:            %s\n"
        "Run Date:             %s\n"
        "Operating System:     %s\n"
        "Input Arguments:\n",
        history->rss_ringoccs_version,
        history->libtmpl_version,
        history->c_version,
        history->host_name,
        history->user_name,
        history->run_date,
        history->operating_system
    );

    if (!history->input_vars[0])
        puts("None");
    else
    {
        const size_t var_len = TMPL_ARRAY_SIZE(history->input_vars);
        size_t n = 0;

        while (history->input_vars[n])
        {
            printf(
                "\t%s: %s\n",
                history->input_vars[n],
                history->input_vars[n+1]
            );

            n += 2;

            if (n >= var_len)
                break;
        }
    }

    puts("Input Keywords:");

    if (!history->input_kwds[0])
        puts("None");
    else
    {
        const size_t var_len = TMPL_ARRAY_SIZE(history->input_kwds);
        size_t n = 0;

        while (history->input_kwds[n])
        {
            printf(
                "\t%s: %s\n",
                history->input_kwds[n],
                history->input_kwds[n+1]
            );

            n += 2;

            if (n >= var_len)
                break;
        }
    }
}
