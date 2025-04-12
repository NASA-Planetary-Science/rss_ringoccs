/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking Tau parameters for possible errors.                 */
void rssringoccs_Tau_Check_Geometry(rssringoccs_TAUObj * const tau)
{
    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    rssringoccs_Tau_Check_Azimuth_Angle(tau);
    rssringoccs_Tau_Check_Displacement(tau);
    rssringoccs_Tau_Check_Opening_Angle(tau);
    rssringoccs_Tau_Check_Ring_Distance(tau);
    rssringoccs_Tau_Check_Ring_Radius(tau);
}
/*  End of rssringoccs_Tau_Check_Geometry.                                    */
