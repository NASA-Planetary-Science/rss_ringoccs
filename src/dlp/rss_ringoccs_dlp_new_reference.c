/*  Header file with the DLPObj typedef.                                      */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  NULL provided here.                                                       */
#include <stddef.h>

/*  Function prototype / forward declaration.                                 */
extern rssringoccs_DLPObj *
rssringoccs_DLP_New_Reference(rssringoccs_DLPObj * const dlp);

/*  Function for safely incrementing the reference counter for a DLP object.  */
rssringoccs_DLPObj *
rssringoccs_DLP_New_Reference(rssringoccs_DLPObj * const dlp)
{
    if (!dlp)
        return NULL;

    if (dlp->error_occurred)
        return NULL;

    ++dlp->reference_count;

    return dlp;
}
/*  End of rssringoccs_DLP_New_Reference.                                     */
