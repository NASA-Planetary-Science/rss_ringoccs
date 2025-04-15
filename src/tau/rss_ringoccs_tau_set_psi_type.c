
/*  String manipulation tools and Booleans provided here.                     */
#include <libtmpl/include/tmpl.h>

/*  Function prototype and Tau object typedef given here.                     */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  The error message for psi-type is longer than the length required by the  *
 *  C89 standard, which is 509 characters. One way to get around this in a    *
 *  portable manner is to use a char array instead of a string literal.       *
 *  The following arrays prints out the following:                            *
 *                                                                            *
 *      Error Encountered: rss_ringoccs                                       *
 *          rssringoccs_Tau_Set_Psi_Type                                      *
 *                                                                            *
 *      Illegal string for psitype. Allowed strings:                          *
 *          newton:     Newton-Raphson method                                 *
 *          newtond:    Newton-Raphson with D perturbation.                   *
 *          newtondold: Newton-Raphson with the old D algorithm.              *
 *          newtondphi: Newton-Raphson with dD/dphi perturbation.             *
 *          simplefft:  A single FFT of the entire data set.                  *
 *          ellipse:    Newton-Raphson with elliptical perturbation.          *
 *          fresnel:    Quadratic Fresnel approximation                       *
 *          fresneln:   Legendre polynomial approximation with 1<n<256        */
static const char rssringoccs_psi_type_error_message[745] = {
    '\n', '\r', 'E', 'r', 'r', 'o', 'r', ' ', 'E', 'n', 'c', 'o', 'u', 'n',
    't', 'e', 'r', 'e', 'd', ':', ' ', 'r', 's', 's', '_', 'r', 'i', 'n', 'g',
    'o', 'c', 'c', 's', '\n', '\r', '\t', 'r', 's', 's', 'r', 'i', 'n', 'g',
    'o', 'c', 'c', 's', '_', 'T', 'a', 'u', '_', 'S', 'e', 't', '_', 'P', 's',
    'i', '_', 'T', 'y', 'p', 'e', '\n', '\n', '\r', 'I', 'l', 'l', 'e', 'g',
    'a', 'l', ' ', 's', 't', 'r', 'i', 'n', 'g', ' ', 'f', 'o', 'r', ' ', 'p',
    's', 'i', 't', 'y', 'p', 'e', '.', ' ', 'A', 'l', 'l', 'o', 'w', 'e', 'd',
    ' ', 's', 't', 'r', 'i', 'n', 'g', 's', ':', '\n', '\r', '\t', 'n', 'e',
    'w', 't', 'o', 'n', ':', ' ', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o',
    'n', '-', 'R', 'a', 'p', 'h', 's', 'o', 'n', ' ', 'm', 'e', 't', 'h', 'o',
    'd', '\n', '\r', '\t', 'n', 'e', 'w', 't', 'o', 'n', 'd', ':', ' ', ' ',
    ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h', 's', 'o',
    'n', ' ', 'w', 'i', 't', 'h', ' ', 'D', ' ', 'p', 'e', 'r', 't', 'u', 'r',
    'b', 'a', 't', 'i', 'o', 'n', '.', '\n', '\r', '\t', 'n', 'e', 'w', 't',
    'o', 'n', 'd', 'o', 'l', 'd', ':', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-',
    'R', 'a', 'p', 'h', 's', 'o', 'n', ' ', 'w', 'i', 't', 'h', ' ', 't', 'h',
    'e', ' ', 'o', 'l', 'd', ' ', 'D', ' ', 'a', 'l', 'g', 'o', 'r', 'i', 't',
    'h', 'm', '.', '\n', '\r', '\t', 'n', 'e', 'w', 't', 'o', 'n', 'd', 'p',
    'h', 'i', ':', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h',
    's', 'o', 'n', ' ', 'w', 'i', 't', 'h', ' ', 'd', 'D', '/', 'd', 'p', 'h',
    'i', ' ', 'p', 'e', 'r', 't', 'u', 'r', 'b', 'a', 't', 'i', 'o', 'n', '.',
    '\n', '\r', '\t', 's', 'i', 'm', 'p', 'l', 'e', 'f', 'f', 't', ':', ' ',
    ' ', 'A', ' ', 's', 'i', 'n', 'g', 'l', 'e', ' ', 'F', 'F', 'T', ' ', 'o',
    'f', ' ', 't', 'h', 'e', ' ', 'e', 'n', 't', 'i', 'r', 'e', ' ', 'd', 'a',
    't', 'a', ' ', 's', 'e', 't', '.', '\n', '\r', '\t', 'q', 'u', 'a', 'd',
    'r', 'a', 't', 'i', 'c', ':', ' ', ' ', 'Q', 'u', 'a', 'd', 'r', 'a', 't',
    'i', 'c', ' ', 'i', 'n', 't', 'e', 'r', 'p', 'o', 'l', 'a', 't', 'i', 'o',
    'n', ' ', 'o', 'f', ' ', 'n', 'e', 'w', 't', 'o', 'n', '-', 'r', 'a', 'p',
    'h', 's', 'o', 'n', '.', '\r', '\t', 'c', 'u', 'b', 'i', 'c', ':', ' ',
    ' ', ' ', ' ', ' ', ' ', 'C', 'u', 'b', 'i', 'c', ' ', 'i', 'n', 't', 'e',
    'r', 'p', 'o', 'l', 'a', 't', 'i', 'o', 'n', ' ', 'o', 'f', ' ', 'n', 'e',
    'w', 't', 'o', 'n', '-', 'r', 'a', 'p', 'h', 's', 'o', 'n', '.', '\r',
    '\t', 'q', 'u', 'a', 'r', 't', 'i', 'c', ':', ' ', ' ', ' ', ' ', 'Q',
    'u', 'a', 'r', 't', 'i', 'c', ' ', 'i', 'n', 't', 'e', 'r', 'p', 'o', 'l',
    'a', 't', 'i', 'o', 'n', ' ', 'o', 'f', ' ', 'n', 'e', 'w', 't', 'o', 'n',
    '-', 'r', 'a', 'p', 'h', 's', 'o', 'n', '.', '\r', '\t', 'q', 'u', 'a',
    'r', 't', 'i', 'c', 'd', ':', ' ', ' ', ' ', 'Q', 'u', 'a', 'r', 't', 'i',
    'c', ' ', 'i', 'n', 't', 'e', 'r', 'p', 'o', 'l', 'a', 't', 'i', 'o', 'n',
    ' ', 'w', 'i', 't', 'h', ' ', 'D', ' ', 'p', 'e', 'r', 't', 'u', 'r', 'b',
    'a', 't', 'i', 'o', 'n', '.', '\r', '\t', 'e', 'l', 'l', 'i', 'p', 's',
    'e', ':', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a',
    'p', 'h', 's', 'o', 'n', ' ', 'w', 'i', 't', 'h', ' ', 'e', 'l', 'l', 'i',
    'p', 't', 'i', 'c', 'a', 'l', ' ', 'p', 'e', 'r', 't', 'u', 'r', 'b', 'a',
    't', 'i', 'o', 'n', '.', '\n', '\r', '\t', 'f', 'r', 'e', 's', 'n', 'e',
    'l', ':', ' ', ' ', ' ', ' ', 'Q', 'u', 'a', 'd', 'r', 'a', 't', 'i', 'c',
    ' ', 'F', 'r', 'e', 's', 'n', 'e', 'l', ' ', 'a', 'p', 'p', 'r', 'o', 'x',
    'i', 'm', 'a', 't', 'i', 'o', 'n', '\n', '\r', '\t', 'f', 'r', 'e', 's',
    'n', 'e', 'l', 'n', ':', ' ', ' ', ' ', 'L', 'e', 'g', 'e', 'n', 'd', 'r',
    'e', ' ', 'p', 'o', 'l', 'y', 'n', 'o', 'm', 'i', 'a', 'l', ' ', 'a', 'p',
    'p', 'r', 'o', 'x', 'i', 'm', 'a', 't', 'i', 'o', 'n', ' ', 'w', 'i', 't',
    'h', ' ', '1', '<', 'n', '<', '2', '5', '6', '\n', '\0'
};

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

    /*  Newton-Raphson method with corrections for the change in the distance *
     *  between the spacecraft and the ring intercept point over the window.  */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtond"))
        tau->psinum = rssringoccs_PsiType_NewtonD;

    /*  Similar to "newtond", but using the old method.                       */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtondold"))
        tau->psinum = rssringoccs_PsiType_NewtonDOld;

    /*  Newton-Raphson method with corrections for dpsi/dphi accounted for.   */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtondphi"))
        tau->psinum = rssringoccs_PsiType_NewtonDPhi;

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

    /*  Quadratic interpolation to the Newton-Raphson method.                 */
    else if (tmpl_String_Are_Equal(tau_psitype, "quadratic"))
        tau->psinum = rssringoccs_PsiType_NewtonQuadratic;

    /*  Quartic interpolation to the Newton-Raphson method.                   */
    else if (tmpl_String_Are_Equal(tau_psitype, "quartic"))
        tau->psinum = rssringoccs_PsiType_NewtonQuartic;

    /*  Quartic interpolation to the Newton-Raphson method with D correction. */
    else if (tmpl_String_Are_Equal(tau_psitype, "quarticd"))
        tau->psinum = rssringoccs_PsiType_NewtonDQuartic;

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
