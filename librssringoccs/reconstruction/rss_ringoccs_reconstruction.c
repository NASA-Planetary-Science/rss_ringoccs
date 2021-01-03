

#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Reconstruction(rssringoccs_DLPObj *dlp,
                                rssringoccs_TAUObj *tau)
{
    rssringoccs_Copy_DLP_Data_To_Tau(dlp, tau);
    rssringoccs_Check_Tau_Data(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Check_Tau_Data_Range(tau);
    rssringoccs_Tau_Select_Window_Func(tau);

    if      (tau->psinum == rssringoccs_DR_Fresnel)
        rssringoccs_Diffraction_Correction_Fresnel(tau);
    else if (tau->psinum == rssringoccs_DR_Legendre)
        rssringoccs_Diffraction_Correction_Legendre(tau);
    else
        rssringoccs_Diffraction_Correction_Newton(tau);
    return;
}
