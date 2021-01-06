#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Tau_Set_WType(const char *wtype, rssringoccs_TAUObj *tau)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (wtype == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_WType\n\n"
            "\rInput string is NULL. Returning.\n"
        );
        return;
    }

    if (tau->wtype != NULL)
        free(tau->wtype);

    tau->wtype = rssringoccs_strdup(wtype);
    rssringoccs_Remove_Spaces(tau->wtype);
    rssringoccs_Make_Lower(tau->wtype);

    if (strcmp(tau->wtype, "rect") == 0)
    {
        tau->normeq = RectNormEQ;
        tau->window_func = rssringoccs_Double_Rect_Window;
    }
    else if (strcmp(tau->wtype, "coss") == 0)
    {
        tau->normeq = CossNormEQ;
        tau->window_func = rssringoccs_Double_Coss_Window;
    }
    else if (strcmp(tau->wtype, "kb20") == 0)
    {
        tau->normeq = KB20NormEQ;
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_2_0;
    }
    else if (strcmp(tau->wtype, "kb25") == 0)
    {
        tau->normeq = KB25NormEQ;
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_2_5;
    }
    else if (strcmp(tau->wtype, "kb35") == 0)
    {
        tau->normeq = KB35NormEQ;
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_3_5;
    }
    else if (strcmp(tau->wtype, "kbmd20") == 0)
    {
        tau->normeq = KBMD20NormEQ;
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_2_0;
    }
    else if (strcmp(tau->wtype, "kbmd25") == 0)
    {
        tau->normeq = KBMD25NormEQ;
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_2_5;
    }
    else if (strcmp(tau->wtype, "kbmd35") == 0)
    {
        tau->normeq = KBMD35NormEQ;
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_3_5;
    }
    else
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_WType\n\n"
            "\r\tIllegal string for wtype. Allowed strings are:\n"
            "\r\t\trect:    Rectangular Window\n"
            "\r\t\tcoss:    Squared Cosine Window\n"
            "\r\t\tkb20:    Kaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkb25:    Kaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkb35:    Kaiser-Bessel with alpha=3.5 pi\n"
            "\r\t\tkbmd20:  Modified Kaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkbmd25:  Modified Kaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkbmd35:  Modified Kaiser-Bessel with alpha=3.5 pi\n"
        );
        tau->normeq = -1.0;
    }
}
