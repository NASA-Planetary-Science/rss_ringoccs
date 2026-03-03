/*  Header file with the DLPObj typedef.                                      */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

void rssringoccs_DLP_Release(rssringoccs_DLPObj *dlp)
{
    if (!dlp)
        return;

    if (dlp->reference_count == 0)
        return;

    --dlp->reference_count;

    if (dlp->reference_count == 0)
        rssringoccs_DLP_Destroy(&dlp);
}
