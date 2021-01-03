#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <stdlib.h>

void rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau)
{
    if (tau == NULL)
        return;
    else if (tau->error_occurred)
        return;

    if (tau->w_km_vals != NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rtau->w_km_vals is not NULL. You have either already processed\n"
            "\rthe window width or did not create the rssringoccs_TAUObj\n"
            "\rpointer with rssringoccs_Create_TAUObj. Returning with error.\n"
        );
        return;
    }
}
