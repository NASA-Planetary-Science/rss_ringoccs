#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Diffraction_Correction(rssringoccs_TAUObj *tau)
{
    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    if (tau->psinum == rssringoccs_PsiType_Fresnel)
        rssringoccs_Diffraction_Correction_Fresnel(tau);

    else if (tau->psinum == rssringoccs_PsiType_Legendre)
        rssringoccs_Diffraction_Correction_Legendre(tau);

    else if (tau->psinum == rssringoccs_PsiType_NewtonSimpleFFT)
        rssringoccs_Diffraction_Correction_SimpleFFT(tau);

    else if (tau->psinum <= rssringoccs_PsiType_NewtonPerturb)
        rssringoccs_Diffraction_Correction_Newton(tau);

    else if (tau->psinum < rssringoccs_PsiType_None)
        rssringoccs_Diffraction_Correction_Newton(tau);

    else
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Reconstruction\n\n"
            "\rpsinum is set to rssringoccs_PsiType_None.\n"
            "\rNo valid psitype requested.\n";
    }
}
