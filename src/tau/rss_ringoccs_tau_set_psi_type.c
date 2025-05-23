
/*  String manipulation tools and Booleans provided here.                     */
#include <libtmpl/include/tmpl.h>

/*  Function prototype and Tau object typedef given here.                     */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

static const char * const rssringoccs_psi_type_error_message =
    "\rError Encountered: rss_ringoccs\n"
    "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
    "\rIllegal string for psitype. Allowed strings:\n"
    "\r\tnewton:     Newton-Raphson method\n"
    "\r\tsimplefft:  A single FFT of the entire data set.\n"
    "\r\tellipse:    Newton-Raphson with elliptical perturbation.\n"
    "\r\tfresnel:    Quadratic Fresnel approximation\n"
    "\r\tfresneln:   Legendre polynomial approximation with 1<n<256\n";

void
rssringoccs_Tau_Set_Psi_Type(const char *psitype, rssringoccs_TAUObj* tau)
{
    char *tau_psitype;

    /*  If the tau pointer is NULL there is nothing to be done.               */
    if (!tau)
        return;

    /*  Similarly if an error occurred before this function was called, abort.*/
    if (tau->error_occurred)
        return;

    /*  If the input string is NULL treat this as an error.                   */
    if (psitype == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
            "\rInput string is NULL. Returning.\n\n";

        return;
    }

    /*  We'll make a copy of the string so that we can remove whitespace and  *
     *  make it lower case without changing the string stored in the Tau      *
     *  object. This allows us to accept case-insensitive inputs.             */
    tau_psitype = tmpl_String_Duplicate(psitype);

    /*  tmpl_String_Duplicate returns NULL on failure. Check for this.        */
    if (!tau_psitype)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
            "\rtmpl_String_Duplicate returned NULL.\n\n";

        return;
    }

    /*  Remove all whitespace from the string. This allows the user to input  *
     *  something like "fresnel 4" instead of "fresnel4" without error.       */
    tmpl_String_Remove_Whitespace(tau_psitype);

    /*  And make all characters lower case. The user can input "Fresnel 4".   */
    tmpl_String_Make_Lower_Case(tau_psitype);

    /*  Standard Newton-Raphson method of finding the stationary azimuth      *
     *  angle for the Fresnel kernel.                                         */
    if (tmpl_String_Are_Equal(tau_psitype, "newton"))
        tau->psinum = rssringoccs_PsiType_Newton;

    /*  Newton-Raphson method, but allowing for arbitrary quartic             *
     *  perturbations to the Fresnel kernel.                                  */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtonperturb"))
        tau->psinum = rssringoccs_PsiType_NewtonPerturb;

    /*  Newton-Rapshon method but for elliptical rings instead of cylindrical.*/
    else if (tmpl_String_Are_Equal(tau_psitype, "ellipse"))
        tau->psinum = rssringoccs_PsiType_NewtonElliptical;

    /*  Computes the reconstruction using a single FFT across the entire data *
     *  set. Only works if the geometry does not very too much across the     *
     *  data. This is possible if the data set is small. Runs in O(N log(N))  *
     *  time, where N is the number of points in the data set, as opposed to  *
     *  the O(N^2) time the usual Newton-Raphson method runs in. Expect very  *
     *  poor results if the range of the data set is very large.              */
    else if (tmpl_String_Are_Equal(tau_psitype, "simplefft"))
        tau->psinum = rssringoccs_PsiType_NewtonSimpleFFT;

    /*  Quartic interpolation to the Newton-Raphson method.                   */
    else if (tmpl_String_Are_Equal(tau_psitype, "quartic"))
        tau->psinum = rssringoccs_PsiType_NewtonQuartic;

    /*  Standard quadratic Fresnel approximation.                             */
    else if (tmpl_String_Are_Equal(tau_psitype, "fresnel"))
        tau->psinum = rssringoccs_PsiType_Fresnel;

    /*  strncmp is a C standard library function that compares the first n    *
     *  elements of two strings. If the first seven elements of tau.psitype   *
     *  are "fresnel", but the string is not exactly "fresnel", try to parse  *
     *  the rest of it and extract a value. For example, if                   *
     *  tau.psitype = "fresnel4", try to extract the "4".                     */
    else if (tmpl_String_Are_First_Characters_Equal(tau_psitype, "fresnel", 7))
    {
        /*  The number part of the string starts after the word "fresnel", so *
         *  start at the 7th character in the string.                         */
        const char *fresnelnum = &tau_psitype[7];

        /*  Convert the string to an integer.                                 */
        tau->order = tmpl_String_To_UChar(fresnelnum);
        tau->psinum = rssringoccs_PsiType_Legendre;

        /*  Zero is not a valid input (like "fresnel0"), and atol returns 0   *
         *  if it could not parse the input. Treat this as an error.          */
        if (tau->order == 0)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
                "\rCould not parse psitype. tmpl_String_To_UChar returned\n"
                "\rzero. Your input has 'fresnel' in it but either has an\n"
                "\rinvalid entry after, or a zero.\n\n";

            tau->order = 0U;
            goto FINISH;
        }

        /*  Fresnel1 is invalid. Check for this.                              */
        else if (tau->order == 1)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
                "\rfresnel1 is not allowed since the constant and linear\n"
                "\rterms of the Fresnel kernel are zero. Please choose\n"
                "\r'fresneln' with n > 1.\n\n";

            tau->order = 0U;
            goto FINISH;
        }

        /*  'fresnel2' is treated as 'fresnel'.                               */
        else if (tau->order == 2)
        {
            tau->psinum = rssringoccs_PsiType_Fresnel;
            tau->order = 0U;
            goto FINISH;
        }

        /*  tau->order stores the number of coefficients needed. The constant *
         *  and linear terms of the expansion are zero, so a degree n         *
         *  expansion needs n - 1 terms. Decrement the order by 1.            */
        --tau->order;
    }

    /*  If we get here we have an invalid string. Set an error.               */
    else
    {
        tau->order = 0;
        tau->psinum = rssringoccs_PsiType_None;
        tau->error_occurred = tmpl_True;
        tau->error_message = rssringoccs_psi_type_error_message;
    }

FINISH:

    /*  tmpl_String_Duplicate calls malloc so we need to free this memory.    */
    tmpl_String_Destroy(&tau_psitype);
}
/*  End of rssringoccs_Tau_Set_Psi_Type.                                      */
