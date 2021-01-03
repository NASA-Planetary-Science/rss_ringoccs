

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    rssringoccs_Check_Tau_Data(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Check_Tau_Data_Range(tau);
    rssringoccs_Tau_Select_Window_Func(tau);

    tau->T_out = malloc(sizeof(*tau->T_out) * tau->arr_size);

    if      (tau->psinum == rssringoccs_DR_Fresnel)
        rssringoccs_Diffraction_Correction_Fresnel(tau);
    else if (tau->psinum == rssringoccs_DR_Legendre)
        rssringoccs_Diffraction_Correction_Legendre(tau);
    else
        rssringoccs_Diffraction_Correction_Newton(tau);
    return;
}
