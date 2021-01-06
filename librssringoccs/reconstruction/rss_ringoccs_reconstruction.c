
#include <stdlib.h>
#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    rssringoccs_Tau_Check_Keywords(tau);
    rssringoccs_Tau_Check_Occ_Type(tau);
    rssringoccs_Tau_Compute_Vars(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Tau_Check_Data_Range(tau);

    tau->T_out = calloc(tau->arr_size, sizeof(*tau->T_out));
    rssringoccs_Tau_Check_Data(tau);

    if      (tau->psinum == rssringoccs_DR_Fresnel)
        rssringoccs_Diffraction_Correction_Fresnel(tau);
    else if (tau->psinum == rssringoccs_DR_Legendre)
        rssringoccs_Diffraction_Correction_Legendre(tau);
    else
        rssringoccs_Diffraction_Correction_Newton(tau);

    rssringoccs_Tau_Finish(tau);
    return;
}
