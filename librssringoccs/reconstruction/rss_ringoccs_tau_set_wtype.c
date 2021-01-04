#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Tau_Set_WType(const char *wtype, rssringoccs_TAUObj *tau)
{
    unsigned long str_len;
    unsigned long new_str_len, n;

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

    str_len = strlen(wtype);
    if (str_len == 0)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_WType\n\n"
            "\rInput string is empty. Returning.\n"
        );
        return;
    }

    tau->wtype = malloc(sizeof(*tau->wtype) * str_len);
    new_str_len = 0;
    n = 0;

    while(wtype[n])
    {
        if (!isspace(wtype[n]))
        {
            tau->wtype[new_str_len] = tolower(wtype[n]);
            ++new_str_len;
        }
        ++n;
    }

    if (new_str_len == 0)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_WType\n\n"
            "\rInput string is only spaces. Returning.\n"
        );
        return;
    }

    tau->wtype = realloc(tau->wtype, new_str_len);

    if (strcmp(tau->wtype, "rect") == 0)
        tau->normeq = RectNormEQ;
    else if (strcmp(tau->wtype, "coss") == 0)
        tau->normeq = CossNormEQ;
    else if (strcmp(tau->wtype, "kb20") == 0)
        tau->normeq = KB20NormEQ;
    else if (strcmp(tau->wtype, "kb25") == 0)
        tau->normeq = KB25NormEQ;
    else if (strcmp(tau->wtype, "kb35") == 0)
        tau->normeq = KB35NormEQ;
    else if (strcmp(tau->wtype, "kbmd20") == 0)
        tau->normeq = KBMD20NormEQ;
    else if (strcmp(tau->wtype, "kbmd25") == 0)
        tau->normeq = KBMD25NormEQ;
    else if (strcmp(tau->wtype, "kbmd35") == 0)
        tau->normeq = KBMD35NormEQ;
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
