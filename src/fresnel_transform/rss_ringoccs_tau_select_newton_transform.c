#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

rssringoccs_FresnelNewtonTransform
rssringoccs_Tau_Select_Newton_Transform(rssringoccs_TAUObj * const tau)
{
    size_t start, num;

    if (!tau)
        return NULL;

    if (tau->error_occurred)
        return NULL;

    if (tau->psinum == rssringoccs_PsiType_None)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Select_Newton_Transform\n\n"
            "\rpsinum is set to None type.\n";

        return NULL;
    }

    if (tau->psinum < rssringoccs_PsiType_NewtonRiemann)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Select_Newton_Transform\n\n"
            "\rpsinum is not a Newton type.\n";

        return NULL;
    }

    if (tau->psinum > rssringoccs_PsiType_NewtonPerturb)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Select_Newton_Transform\n\n"
            "\rpsinum is set to an interpolated Newton transform,\n"
            "\rnot a Newton transform.\n";

        return NULL;
    }

    num = TMPL_CAST(tau->psinum, size_t);
    start = TMPL_CAST(rssringoccs_PsiType_NewtonRiemann, size_t);
    num -= start;

    return rssringoccs_newton_transform_table[num];
}
