/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Function for checking DLP parameters for possible errors.                 */
void rssringoccs_DLP_Check_Geometry(rssringoccs_DLPObj * const dlp)
{
    if (!dlp)
        return;

    if (dlp->error_occurred)
        return;

    rssringoccs_DLP_Check_Azimuth_Angle(dlp);
    rssringoccs_DLP_Check_Displacement(dlp);
    rssringoccs_DLP_Check_Opening_Angle(dlp);
    rssringoccs_DLP_Check_Ring_Distance(dlp);
    rssringoccs_DLP_Check_Ring_Radius(dlp);
}
/*  End of rssringoccs_DLP_Check_Geometry.                                    */
