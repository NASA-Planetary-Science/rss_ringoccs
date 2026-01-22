/*  TMPL_RESTRICT macro provided here.                                        */
#include <libtmpl/include/tmpl_config.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  String manipulation functions found here.                                 */
#include <libtmpl/include/tmpl_string.h>

/*  The psitype enum, with all possible psitypes, is found here.              */
#include <rss_ringoccs/include/types/rss_ringoccs_psitype.h>

/*  Tau object typedef given here.                                            */
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>

/*  Function prototype / forward declaration.                                 */
extern void
rssringoccs_Tau_Set_Psi_Type(const char * TMPL_RESTRICT const psitype,
                             rssringoccs_TAUObj * TMPL_RESTRICT const tau);

/*  Error message for the legal psitypes. Defined at the bottom of this file. */
static const char rssringoccs_psi_type_error_message[2708];

/*  Function for setting the psitype in a Tau object from a string.           */
void
rssringoccs_Tau_Set_Psi_Type(const char * TMPL_RESTRICT const psitype,
                             rssringoccs_TAUObj * TMPL_RESTRICT const tau)
{
    char *tau_psitype;

    /*  If the tau pointer is NULL there is nothing to be done.               */
    if (!tau)
        return;

    /*  Similarly if an error occurred before this function was called, abort.*/
    if (tau->error_occurred)
        return;

    /*  If the input string is NULL treat this as an error.                   */
    if (!psitype)
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
    if (tmpl_String_Are_Equal(tau_psitype, "newtonriemann"))
        tau->psinum = rssringoccs_PsiType_NewtonRiemann;

    else if (tmpl_String_Are_Equal(tau_psitype, "newtonfilon01"))
        tau->psinum = rssringoccs_PsiType_NewtonFilon01;

    else if (tmpl_String_Are_Equal(tau_psitype, "newtonfilon11"))
        tau->psinum = rssringoccs_PsiType_NewtonFilon11;

    else if (tmpl_String_Are_Equal(tau_psitype, "newtonfilon02"))
        tau->psinum = rssringoccs_PsiType_NewtonFilon02;

    else if (tmpl_String_Are_Equal(tau_psitype, "newtonfilon12"))
        tau->psinum = rssringoccs_PsiType_NewtonFilon12;

    /*  Newton-Raphson method, but allowing for arbitrary quartic             *
     *  perturbations to the Fresnel kernel.                                  */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtonperturb"))
        tau->psinum = rssringoccs_PsiType_NewtonPerturb;

    /*  Newton-Rapshon method but for elliptical rings instead of cylindrical.*/
    else if (tmpl_String_Are_Equal(tau_psitype, "ellipse"))
        tau->psinum = rssringoccs_PsiType_EllipticNewton;

    /*  Computes the reconstruction using a single FFT across the entire data *
     *  set. Only works if the geometry does not very too much across the     *
     *  data. This is possible if the data set is small. Runs in O(N log(N))  *
     *  time, where N is the number of points in the data set, as opposed to  *
     *  the O(N^2) time the usual Newton-Raphson method runs in. Expect very  *
     *  poor results if the range of the data set is very large.              */
    else if (tmpl_String_Are_Equal(tau_psitype, "simplefft"))
        tau->psinum = rssringoccs_PsiType_NewtonSimpleFFT;

    /*  Quartic interpolation to the Newton-Raphson method.                   */
    else if (tmpl_String_Are_Equal(tau_psitype, "newton4"))
        tau->psinum = rssringoccs_PsiType_Newton4;

    /*  Octic interpolation to the Newton-Raphson method.                     */
    else if (tmpl_String_Are_Equal(tau_psitype, "newton8"))
        tau->psinum = rssringoccs_PsiType_Newton8;

    /*  Octic interpolation to the Newton-Raphson method.                     */
    else if (tmpl_String_Are_Equal(tau_psitype, "newton16"))
        tau->psinum = rssringoccs_PsiType_Newton16;

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

            goto FINISH;
        }

        /*  Fresnel1 is invalid. Check for this.                              */
        if (tau->order == 1)
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
        if (tau->order == 2)
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

/*  The string for the error message. C89 compilers are required to support   *
 *  string literals of up to ~500 characters. The psitype error message is    *
 *  over 2000 characters long, so we provide it as a char array. This makes   *
 *  it infinitely less readable, but portable. The array below was created    *
 *  using a short C program, it was not typed out by hand.                    */
static const char rssringoccs_psi_type_error_message[2708] = {
    '\n', 'E', 'r', 'r', 'o', 'r', ' ', 'E', 'n', 'c', 'o', 'u', 'n', 't',
    'e', 'r', 'e', 'd', ':', ' ', 'r', 's', 's', '_', 'r', 'i', 'n', 'g', 'o',
    'c', 'c', 's', '\n', ' ', ' ', ' ', ' ', 'r', 's', 's', 'r', 'i', 'n',
    'g', 'o', 'c', 'c', 's', '_', 'T', 'a', 'u', '_', 'S', 'e', 't', '_', 'P',
    's', 'i', '_', 'T', 'y', 'p', 'e', '\n', '\n', 'I', 'l', 'l', 'e', 'g',
    'a', 'l', ' ', 's', 't', 'r', 'i', 'n', 'g', ' ', 'f', 'o', 'r', ' ', 'p',
    's', 'i', 't', 'y', 'p', 'e', '.', ' ', 'A', 'l', 'l', 'o', 'w', 'e', 'd',
    ' ', 's', 't', 'r', 'i', 'n', 'g', 's', ':', '\n', ' ', ' ', ' ', ' ',
    'N', 'e', 'w', 't', 'o', 'n', ' ', 'R', 'i', 'e', 'm', 'a', 'n', 'n', ':',
    '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o',
    'n', '-', 'R', 'a', 'p', 'h', 's', 'o', 'n', ' ', 'm', 'e', 't', 'h', 'o',
    'd', ' ', 'u', 's', 'i', 'n', 'g', ' ', 'a', ' ', 'R', 'i', 'e', 'm', 'a',
    'n', 'n', ' ', 's', 'u', 'm', ' ', 'f', 'o', 'r', ' ', 't', 'h', 'e', ' ',
    'i', 'n', 't', 'e', 'g', 'r', 'a', 'l', '.', '\n', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', 'T', 'h', 'i', 's', ' ', 'i', 's', ' ', 'i', 'n', 'a',
    'c', 'c', 'u', 'r', 'a', 't', 'e', ' ', 'f', 'o', 'r', ' ', 'd', 'a', 't',
    'a', ' ', 's', 'e', 't', 's', ' ', 'w', 'h', 'e', 'r', 'e', ' ', 'd', ' ',
    'p', 's', 'i', ' ', '/', ' ', 'd', ' ', 'r', 'h', 'o', ' ', 'i', 's', ' ',
    'l', 'a', 'r', 'g', 'e', ',', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', 'a', 'n', 'd', ' ', 't', 'h', 'e', ' ', 'F', 'i', 'l', 'o', 'n', ' ',
    'm', 'e', 't', 'h', 'o', 'd', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b',
    'e', ' ', 'u', 's', 'e', 'd', ' ', 'i', 'n', 's', 't', 'e', 'a', 'd', '.',
    '\n', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ', 'F', 'i',
    'l', 'o', 'n', ' ', '0', '1', ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h', 's', 'o',
    'n', ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'u', 's', 'i', 'n', 'g', ' ',
    'a', ' ', 'm', 'o', 'd', 'i', 'f', 'i', 'e', 'd', ' ', 'l', 'i', 'n', 'e',
    'a', 'r', ' ', 'F', 'i', 'l', 'o', 'n', ' ', 'q', 'u', 'a', 'd', 'r', 'a',
    't', 'u', 'r', 'e', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'f',
    'o', 'r', ' ', 't', 'h', 'e', ' ', 'i', 'n', 't', 'e', 'g', 'r', 'a', 'l',
    '.', ' ', 'T', 'h', 'e', ' ', 'd', 'a', 't', 'a', ' ', '\'', 'T', '\'',
    ' ', 'i', 's', ' ', 'a', 's', 's', 'u', 'm', 'e', 'd', ' ', 'c', 'o', 'n',
    's', 't', 'a', 'n', 't', ' ', 'a', 'c', 'r', 'o', 's', 's', ' ', 'e', 'a',
    'c', 'h', ' ', 'b', 'i', 'n', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', 'a', 'n', 'd', ' ', 't', 'h', 'e', ' ', 'F', 'r', 'e', 's', 'n', 'e',
    'l', ' ', 'p', 'h', 'a', 's', 'e', ' ', '\'', 'p', 's', 'i', '\'', ' ',
    'i', 's', ' ', 'f', 'i', 't', 't', 'e', 'd', ' ', 't', 'o', ' ', 'a', ' ',
    'l', 'i', 'n', 'e', ' ', 'p', 's', 'i', ' ', '=', ' ', 'a', 'x', ' ', '+',
    ' ', 'b', '.', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T', 'h',
    'a', 't', ' ', 'i', 's', ',', ' ', 'T', ' ', 'u', 's', 'e', 's', ' ', 'a',
    ' ', 'd', 'e', 'g', 'r', 'e', 'e', ' ', '0', ' ', 'a', 'p', 'p', 'r', 'o',
    'x', 'i', 'm', 'a', 't', 'i', 'o', 'n', ' ', 'a', 'n', 'd', ' ', 'p', 's',
    'i', ' ', 'u', 's', 'e', 's', ' ', 'd', 'e', 'g', 'r', 'e', 'e', ' ', '1',
    '.', '\n', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ', 'F',
    'i', 'l', 'o', 'n', ' ', '1', '1', ':', '\n', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h', 's',
    'o', 'n', ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'u', 's', 'i', 'n', 'g',
    ' ', 'a', ' ', 'm', 'o', 'd', 'i', 'f', 'i', 'e', 'd', ' ', 'l', 'i', 'n',
    'e', 'a', 'r', ' ', 'F', 'i', 'l', 'o', 'n', ' ', 'q', 'u', 'a', 'd', 'r',
    'a', 't', 'u', 'r', 'e', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    'f', 'o', 'r', ' ', 't', 'h', 'e', ' ', 'i', 'n', 't', 'e', 'g', 'r', 'a',
    'l', '.', ' ', 'T', 'h', 'e', ' ', 'd', 'a', 't', 'a', ' ', '\'', 'T',
    '\'', ' ', 'a', 'n', 'd', ' ', 't', 'h', 'e', ' ', 'F', 'r', 'e', 's',
    'n', 'e', 'l', ' ', 'p', 'h', 'a', 's', 'e', ' ', '\'', 'p', 's', 'i',
    '\'', ' ', 'a', 'r', 'e', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    'b', 'o', 't', 'h', ' ', 'a', 's', 's', 'u', 'm', 'e', 'd', ' ', 'l', 'i',
    'n', 'e', 'a', 'r', ' ', 'a', 'c', 'r', 'o', 's', 's', ' ', 'a', ' ', 'b',
    'i', 'n', ',', ' ', 'T', ' ', '=', ' ', 'a', 'x', ' ', '+', ' ', 'b', ',',
    ' ', 'p', 's', 'i', ' ', '=', ' ', 'c', 'x', ' ', '+', ' ', 'd', '.', '\n',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T', 'h', 'a', 't', ' ', 'i', 's',
    ',', ' ', 'T', ' ', 'u', 's', 'e', 's', ' ', 'a', ' ', 'd', 'e', 'g', 'r',
    'e', 'e', ' ', '1', ' ', 'a', 'p', 'p', 'r', 'o', 'x', 'i', 'm', 'a', 't',
    'i', 'o', 'n', ' ', 'a', 'n', 'd', ' ', 'p', 's', 'i', ' ', 'u', 's', 'e',
    's', ' ', 'd', 'e', 'g', 'r', 'e', 'e', ' ', '1', '.', '\n', ' ', ' ',
    ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ', 'P', 'e', 'r', 't', 'u', 'r',
    'b', ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'S', 'a', 'm',
    'e', ' ', 'a', 's', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ', 'Q', 'u', 'a',
    'd', 'r', 'a', 't', 'i', 'c', ' ', 'F', 'i', 'l', 'o', 'n', ',', ' ', 'b',
    'u', 't', ' ', 'a', 'n', ' ', 'a', 'r', 'b', 'i', 't', 'r', 'a', 'r', 'y',
    ' ', 'q', 'u', 'a', 'r', 't', 'i', 'c', ' ', 'i', 's', '\n', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', 'a', 'l', 'l', 'o', 'w', 'e', 'd', ' ', 'f',
    'o', 'r', ' ', 'p', 'e', 'r', 't', 'u', 'r', 'b', 'i', 'n', 'g', ' ', 't',
    'h', 'e', ' ', 's', 't', 'a', 't', 'i', 'o', 'n', 'a', 'r', 'y', ' ', 'F',
    'r', 'e', 's', 'n', 'e', 'l', ' ', 'p', 'h', 'a', 's', 'e', ' ', '(', 'p',
    's', 'i', ')', '.', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T',
    'h', 'i', 's', ' ', 'm', 'a', 'y', ' ', 'b', 'e', ' ', 's', 'e', 't', ' ',
    'a', 't', ' ', 't', 'h', 'e', ' ', 'p', 'y', 't', 'h', 'o', 'n', ' ', 'l',
    'e', 'v', 'e', 'l', ' ', 'u', 's', 'i', 'n', 'g', ' ', 't', 'h', 'e', ' ',
    'p', 'e', 'r', 't', 'u', 'r', 'b', ' ', 'k', 'e', 'y', 'w', 'o', 'r', 'd',
    '.', '\n', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ', '4',
    ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't',
    'o', 'n', '-', 'R', 'a', 'p', 'h', 's', 'o', 'n', ' ', 'm', 'e', 't', 'h',
    'o', 'd', ' ', 'u', 's', 'i', 'n', 'g', ' ', 'a', ' ', 'q', 'u', 'a', 'd',
    'r', 'a', 't', 'i', 'c', ' ', 'F', 'i', 'l', 'o', 'n', ' ', 'q', 'u', 'a',
    'd', 'r', 'a', 't', 'u', 'r', 'e', ' ', 'w', 'i', 't', 'h', ' ', 'a', '\n',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'q', 'u', 'a', 'r', 't', 'i', 'c',
    ' ', 'i', 'n', 't', 'e', 'r', 'p', 'o', 'l', 'a', 't', 'i', 'o', 'n', ' ',
    'o', 'f', ' ', 't', 'h', 'e', ' ', 'F', 'r', 'e', 's', 'n', 'e', 'l', ' ',
    'p', 'h', 'a', 's', 'e', ' ', 'a', 'c', 'r', 'o', 's', 's', ' ', 't', 'h',
    'e', ' ', 'w', 'i', 'n', 'd', 'o', 'w', '.', '\n', ' ', ' ', ' ', ' ',
    'N', 'e', 'w', 't', 'o', 'n', ' ', '8', ':', '\n', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h',
    's', 'o', 'n', ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'u', 's', 'i', 'n',
    'g', ' ', 'a', ' ', 'q', 'u', 'a', 'd', 'r', 'a', 't', 'i', 'c', ' ', 'F',
    'i', 'l', 'o', 'n', ' ', 'q', 'u', 'a', 'd', 'r', 'a', 't', 'u', 'r', 'e',
    ' ', 'w', 'i', 't', 'h', ' ', 'a', 'n', '\n', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', 'o', 'c', 't', 'i', 'c', ' ', 'i', 'n', 't', 'e', 'r', 'p',
    'o', 'l', 'a', 't', 'i', 'o', 'n', ' ', 'o', 'f', ' ', 't', 'h', 'e', ' ',
    'F', 'r', 'e', 's', 'n', 'e', 'l', ' ', 'p', 'h', 'a', 's', 'e', ' ', 'a',
    'c', 'r', 'o', 's', 's', ' ', 't', 'h', 'e', ' ', 'w', 'i', 'n', 'd', 'o',
    'w', '.', '\n', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', ' ',
    '1', '6', ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'N', 'e',
    'w', 't', 'o', 'n', '-', 'R', 'a', 'p', 'h', 's', 'o', 'n', ' ', 'm', 'e',
    't', 'h', 'o', 'd', ' ', 'u', 's', 'i', 'n', 'g', ' ', 'a', ' ', 'q', 'u',
    'a', 'd', 'r', 'a', 't', 'i', 'c', ' ', 'F', 'i', 'l', 'o', 'n', ' ', 'q',
    'u', 'a', 'd', 'r', 'a', 't', 'u', 'r', 'e', ' ', 'w', 'i', 't', 'h', ' ',
    'a', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'h', 'e', 'x', 'a',
    'd', 'e', 'c', 'i', 'c', ' ', 'i', 'n', 't', 'e', 'r', 'p', 'o', 'l', 'a',
    't', 'i', 'o', 'n', ' ', 'o', 'f', ' ', 't', 'h', 'e', ' ', 'F', 'r', 'e',
    's', 'n', 'e', 'l', ' ', 'p', 'h', 'a', 's', 'e', ' ', 'a', 'c', 'r', 'o',
    's', 's', ' ', 't', 'h', 'e', ' ', 'w', 'i', 'n', 'd', 'o', 'w', '.', '\n',
    ' ', ' ', ' ', ' ', 'S', 'i', 'm', 'p', 'l', 'e', ' ', 'F', 'F', 'T', ':',
    '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'U', 's', 'e', ' ', 'a',
    ' ', 's', 'i', 'n', 'g', 'l', 'e', ' ', 'F', 'F', 'T', ' ', 'o', 'f', ' ',
    't', 'h', 'e', ' ', 'e', 'n', 't', 'i', 'r', 'e', ' ', 'd', 'a', 't', 'a',
    ' ', 's', 'e', 't', '.', ' ', 'T', 'h', 'i', 's', ' ', 'i', 's', ' ', 'v',
    'e', 'r', 'y', ' ', 'f', 'a', 's', 't', ',', ' ', 'b', 'u', 't', '\n',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'i', 't', ' ', 'a', 's', 's', 'u',
    'm', 'e', 's', ' ', 't', 'h', 'e', ' ', 'g', 'e', 'o', 'm', 'e', 't', 'r',
    'y', ' ', 'i', 's', ' ', 'n', 'e', 'a', 'r', 'l', 'y', ' ', 'c', 'o', 'n',
    's', 't', 'a', 'n', 't', ' ', 'f', 'o', 'r', ' ', 't', 'h', 'e', ' ', 'e',
    'n', 't', 'i', 'r', 'e', ' ', 'd', 'a', 't', 'a', '\n', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', 's', 'e', 't', '.', ' ', 'O', 'n', 'l', 'y', ' ',
    'u', 's', 'e', ' ', 't', 'h', 'i', 's', ' ', 'f', 'o', 'r', ' ', 's', 'm',
    'a', 'l', 'l', ' ', 'r', 'a', 'n', 'g', 'e', 's', ',', ' ', 'o', 'r', ' ',
    's', 'e', 't', 's', ' ', 'w', 'i', 't', 'h', ' ', 's', 't', 'a', 'b', 'l',
    'e', ' ', 'g', 'e', 'o', 'm', 'e', 't', 'r', 'y', '.', '\n', ' ', ' ',
    ' ', ' ', 'E', 'l', 'l', 'i', 'p', 's', 'e', ':', '\n', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', 'N', 'e', 'w', 't', 'o', 'n', '-', 'R', 'a', 'p',
    'h', 's', 'o', 'n', ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'u', 's', 'i',
    'n', 'g', ' ', 't', 'h', 'e', ' ', 'e', 'l', 'l', 'i', 'p', 't', 'i', 'c',
    'a', 'l', ' ', 'p', 'e', 'r', 't', 'u', 'r', 'b', 'a', 't', 'i', 'o', 'n',
    '.', ' ', 'T', 'h', 'i', 's', ' ', 'i', 's', '\n', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', 'u', 's', 'e', 'f', 'u', 'l', ' ', 'f', 'o', 'r', ' ',
    'e', 'l', 'l', 'i', 'p', 't', 'i', 'c', 'a', 'l', ' ', 'r', 'i', 'n', 'g',
    's', '.', '\n', ' ', ' ', ' ', ' ', 'F', 'r', 'e', 's', 'n', 'e', 'l',
    ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'C', 'l', 'a', 's',
    's', 'i', 'c', ' ', 'q', 'u', 'a', 'd', 'r', 'a', 't', 'i', 'c', ' ', 'F',
    'r', 'e', 's', 'n', 'e', 'l', ' ', 'a', 'p', 'p', 'r', 'o', 'x', 'i', 'm',
    'a', 't', 'i', 'o', 'n', '.', ' ', 'T', 'h', 'i', 's', ' ', 'i', 's', ' ',
    't', 'h', 'e', ' ', 'f', 'a', 's', 't', 'e', 's', 't', '\n', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', 'm', 'e', 't', 'h', 'o', 'd', ',', ' ', 'b',
    'u', 't', ' ', 'i', 't', ' ', 'i', 's', ' ', 'n', 'o', 't', ' ', 'a', 'p',
    'p', 'l', 'i', 'c', 'a', 'b', 'l', 'e', ' ', 't', 'o', ' ', 'm', 'a', 'n',
    'y', ' ', 'd', 'a', 't', 'a', ' ', 's', 'e', 't', 's', '.', '\n', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'R', 'e', 'v', '0', '0', '7', ' ', 'f',
    'r', 'o', 'm', ' ', 'C', 'a', 's', 's', 'i', 'n', 'i', ' ', 'a', 'n', 'd',
    ' ', 't', 'h', 'e', ' ', 'V', 'o', 'y', 'a', 'g', 'e', 'r', ' ', '2', ' ',
    'U', 'r', 'a', 'n', 'u', 's', ' ', 'd', 'a', 't', 'a', ' ', 'a', 'r', 'e',
    ' ', 'e', 'x', 'a', 'm', 'p', 'l', 'e', 's', '\n', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', 'o', 'f', ' ', 'i', 'd', 'e', 'a', 'l', ' ', 's', 'e',
    't', 's', '.', ' ', 'F', 'o', 'r', ' ', 'C', 'a', 's', 's', 'i', 'n', 'i',
    '\'', 's', ' ', 'R', 'e', 'v', '1', '3', '3', ' ', 't', 'h', 'i', 's',
    ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'i', 's', ' ', 'u', 'n', 'u', 's',
    'a', 'b', 'l', 'e', '.', '\n', ' ', ' ', ' ', ' ', 'F', 'r', 'e', 's',
    'n', 'e', 'l', ' ', 'n', ':', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', 'L', 'e', 'g', 'e', 'n', 'd', 'r', 'e', ' ', 'p', 'o', 'l', 'y', 'n',
    'o', 'm', 'i', 'a', 'l', ' ', 'a', 'p', 'p', 'r', 'o', 'x', 'i', 'm', 'a',
    't', 'i', 'o', 'n', ' ', 'w', 'i', 't', 'h', ' ', '3', ' ', '<', ' ', 'n',
    ' ', '<', ' ', '2', '5', '6', '.', ' ', 'L', 'i', 'k', 'e', ' ', 't', 'h',
    'e', '\n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'F', 'r', 'e', 's',
    'n', 'e', 'l', ' ', 'm', 'e', 't', 'h', 'o', 'd', ',', ' ', 't', 'h', 'i',
    's', ' ', 'i', 's', ' ', 'o', 'n', 'l', 'y', ' ', 'a', 'p', 'p', 'l', 'i',
    'c', 'a', 'b', 'l', 'e', ' ', 't', 'o', ' ', 'c', 'e', 'r', 't', 'a', 'i',
    'n', ' ', 'd', 'a', 't', 'a', ' ', 's', 'e', 't', 's', '.', '\n', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'U', 's', 'e', ' ', 't', 'h', 'i', 's',
    ' ', 'w', 'h', 'e', 'n', ' ', 't', 'h', 'e', ' ', 'o', 'p', 'e', 'n', 'i',
    'n', 'g', ' ', 'a', 'n', 'g', 'l', 'e', ' ', 'i', 's', ' ', 'l', 'a', 'r',
    'g', 'e', '.', ' ', 'I', 't', ' ', 'i', 's', ' ', 'n', 'e', 'a', 'r', 'l',
    'y', ' ', 'a', 's', ' ', 'f', 'a', 's', 't', ' ', 'a', 's', '\n', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', 't', 'h', 'e', ' ', 'F', 'r', 'e', 's',
    'n', 'e', 'l', ' ', 'm', 'e', 't', 'h', 'o', 'd', ',', ' ', 'b', 'u', 't',
    ' ', 'm', 'o', 'r', 'e', ' ', 'a', 'c', 'c', 'u', 'r', 'a', 't', 'e', ' ',
    'w', 'h', 'e', 'n', ' ', 'i', 't', ' ', 'i', 's', ' ', 'a', 'p', 'p', 'l',
    'i', 'c', 'a', 'b', 'l', 'e', '.', '\n', '\n', 'N', 'o', 't', 'e', ':',
    '\n', ' ', ' ', ' ', ' ', 'I', 'n', 'p', 'u', 't', ' ', 's', 't', 'r',
    'i', 'n', 'g', 's', ' ', 'a', 'r', 'e', ' ', 'n', 'e', 'i', 't', 'h', 'e',
    'r', ' ', 'c', 'a', 's', 'e', ' ', 'n', 'o', 'r', ' ', 's', 'p', 'a', 'c',
    'e', ' ', 's', 'e', 'n', 's', 'i', 't', 'i', 'v', 'e', '.', ' ', 'Y', 'o',
    'u', ' ', 'm', 'a', 'y', ' ', 'a', 'd', 'd', '\n', ' ', ' ', ' ', ' ',
    'o', 'r', ' ', 'r', 'e', 'm', 'o', 'v', 'e', ' ', 's', 'p', 'a', 'c', 'e',
    's', ',', ' ', 'a', 'n', 'd', ' ', 'u', 's', 'e', ' ', 'u', 'p', 'p', 'e',
    'r', ' ', 'c', 'a', 's', 'e', ' ', 'o', 'r', ' ', 'l', 'o', 'w', 'e', 'r',
    ' ', 'c', 'a', 's', 'e', ' ', 'l', 'e', 't', 't', 'e', 'r', 's', '.', '\n',
    '\0'
};
