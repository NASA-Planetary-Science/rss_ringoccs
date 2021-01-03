#include <stdlib.h>
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  This function acts like a Python dictionary and returns the correct       *
 *  normalized equivalent width for the corresponding input string. If an     *
 *  invalid string is given this returns -1. The user of this function then   *
 *  has the responsibility of handling this value appropriately, no errors    *
 *  are raised directly from this routine.                                    *
 *                                                                            *
 *  A small note on strcmp (string-compare). This does not simply return 0 or *
 *  1 depending on whether the two inputs are identical, but rather can       *
 *  return negative, positive, or zero, depending on whether one string is    *
 *  "greater" than the other. For identical strings this returns 0.           *
 *  Hence DO NOT write something like:                                        *
 *      if (strcmp(string, same_string))                                      *
 *          do stuff;                                                         *
 *  since strcmp(string, same_string) will return 0, meaning the if-statement *
 *  is false. Instead, do the following:                                      *
 *      if (strcmp(string, same_string) == 0)                                 *
 *          do stuff;                                                         */
void
rssringoccs_Tau_Get_Normeq_From_String(const char *wtype,
                                       rssringoccs_TAUObj *tau)
{
    if (tau == NULL)
        return;

    if (tau->wtype != NULL)
        free(tau->wtype);

    tau->wtype = rssringoccs_strdup(wtype);

    if (strcmp(wtype, "rect") == 0)
        tau->normeq = RectNormEQ;
    else if (strcmp(wtype, "coss") == 0)
        tau->normeq = CossNormEQ;
    else if (strcmp(wtype, "kb20") == 0)
        tau->normeq = KB20NormEQ;
    else if (strcmp(wtype, "kb25") == 0)
        tau->normeq = KB25NormEQ;
    else if (strcmp(wtype, "kb35") == 0)
        tau->normeq = KB35NormEQ;
    else if (strcmp(wtype, "kbmd20") == 0)
        tau->normeq = KBMD20NormEQ;
    else if (strcmp(wtype, "kbmd25") == 0)
        tau->normeq = KBMD25NormEQ;
    else if (strcmp(wtype, "kbmd35") == 0)
        tau->normeq = KBMD35NormEQ;
    else
        tau->normeq = -1.0;
}
