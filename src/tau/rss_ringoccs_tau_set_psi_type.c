
/*  String manipulation tools and Booleans provided here.                     */
#include <libtmpl/include/tmpl.h>

/*  Function prototype and Tau object typedef given here.                     */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Error message for an invalid input string.                                */
static const char rssringoccs_psi_err_mes[744] =
    "\n\rError Encountered: rss_ringoccs\n"
    "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
    "\rIllegal string for psitype. Allowed strings:\n"
    "\r\tnewton:     Newton-Raphson method\n"
    "\r\tnewtond:    Newton-Raphson with D perturbation.\n"
    "\r\tnewtondold: Newton-Raphson with the old D algorithm.\n"
    "\r\tnewtondphi: Newton-Raphson with dD/dphi perturbation.\n"
    "\r\tsimplefft:  A single FFT of the entire data set.\n"
    "\r\tquadratic:  Quadratic interpolation of newton-raphson."
    "\r\tcubic:      Cubic interpolation of newton-raphson."
    "\r\tquartic:    Quartic interpolation of newton-raphson."
    "\r\tquarticd:   Quartic interpolation with D perturbation."
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
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
            "\rInput string is NULL. Returning.\n\n"
        );

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
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
            "\rtmpl_String_Duplicate returned NULL. Returning.\n\n"
        );

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
        tau->psinum = rssringoccs_DR_Newton;

    /*  Newton-Raphson method with corrections for the change in the distance *
     *  between the spacecraft and the ring intercept point over the window.  */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtond"))
        tau->psinum = rssringoccs_DR_NewtonD;

    /*  Similar to "newtond", but using the old method.                       */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtondold"))
        tau->psinum = rssringoccs_DR_NewtonDOld;

    /*  Newton-Raphson method with corrections for dpsi/dphi accounted for.   */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtondphi"))
        tau->psinum = rssringoccs_DR_NewtonDPhi;

    /*  Newton-Raphson method, but allowing for arbitrary quartic             *
     *  perturbations to the Fresnel kernel.                                  */
    else if (tmpl_String_Are_Equal(tau_psitype, "newtonperturb"))
        tau->psinum = rssringoccs_DR_NewtonPerturb;

    /*  Newton-Rapshon method but for elliptical rings instead of cylindrical.*/
    else if (tmpl_String_Are_Equal(tau_psitype, "ellipse"))
        tau->psinum = rssringoccs_DR_NewtonElliptical;

    /*  Computes the reconstruction using a single FFT across the entire data *
     *  set. Only works if the geometry does not very too much across the     *
     *  data. This is possible if the data set is small. Runs in O(N log(N))  *
     *  time, where N is the number of points in the data set, as opposed to  *
     *  the O(N^2) time the usual Newton-Raphson method runs in. Expect very  *
     *  poor results if the range of the data set is very large.              */
    else if (tmpl_String_Are_Equal(tau_psitype, "simplefft"))
        tau->psinum = rssringoccs_DR_NewtonSimpleFFT;

    /*  Quadratic interpolation to the Newton-Raphson method.                 */
    else if (tmpl_String_Are_Equal(tau_psitype, "quadratic"))
        tau->psinum = rssringoccs_DR_NewtonQuadratic;

    /*  Quartic interpolation to the Newton-Raphson method.                   */
    else if (tmpl_String_Are_Equal(tau_psitype, "quartic"))
        tau->psinum = rssringoccs_DR_NewtonQuartic;

    /*  Quartic interpolation to the Newton-Raphson method with D correction. */
    else if (tmpl_String_Are_Equal(tau_psitype, "quarticd"))
        tau->psinum = rssringoccs_DR_NewtonDQuartic;

    /*  Standard quadratic Fresnel approximation.                             */
    else if (tmpl_String_Are_Equal(tau_psitype, "fresnel"))
        tau->psinum = rssringoccs_DR_Fresnel;

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
        tau->psinum = rssringoccs_DR_Legendre;

        /*  Zero is not a valid input (like "fresnel0"), and atol returns 0   *
         *  if it could not parse the input. Treat this as an error.          */
        if (tau->order == 0)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_String_Duplicate(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Set_Psi_Type\n\n"
                "\rCould not parse psitype. tmpl_String_To_UChar returned\n"
                "\rzero. Your input has 'fresnel' in it but either has an\n"
                "\rinvalid entry after, or a zero.\n\n"
            );

            return;
        }
    }

    /*  If we get here we have an invalid string. Set an error.               */
    else
    {
        tau->order = 0;
        tau->psinum = rssringoccs_DR_None;
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(rssringoccs_psi_err_mes);
    }

    /*  tmpl_String_Duplicate calls malloc so we need to free this memory.    */
    tmpl_String_Destroy(&tau_psitype);
}
/*  End of rssringoccs_Tau_Set_Psi_Type.                                      */
