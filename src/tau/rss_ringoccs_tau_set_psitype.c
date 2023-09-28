#include <stdlib.h>

/*  To help manipulate strings, both native C strings const char* and char*,  *
 *  as well as Python strings passed by the user, include string.h. This is a *
 *  standard library header containing strcmp (string-compare), strncmp       *
 *  (compares first n characters of string), strcpy (string-copy), strncpy,   *
 *  and strlen (string-length), all of which are used at some point.          */
#include <string.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

void
rssringoccs_Tau_Set_Psitype(const char *psitype, rssringoccs_TAUObj* tau)
{
    char *tau_psitype;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (psitype == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psitype\n\n"
            "\rInput string is NULL. Returning.\n"
        );
        return;
    }

    tau_psitype = tmpl_String_Duplicate(psitype);

    if (tau_psitype == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psitype\n\n"
            "\rtmpl_String_Duplicate returned NULL. Returning.\n"
        );
        return;
    }

    tmpl_String_Remove_Whitespace(tau_psitype);
    tmpl_String_Make_Lower_Case(tau_psitype);

    if (strcmp(tau_psitype, "newton") == 0)
        tau->psinum = rssringoccs_DR_Newton;
    else if (strcmp(tau_psitype, "newtond") == 0)
        tau->psinum = rssringoccs_DR_NewtonD;
    else if (strcmp(tau_psitype, "newtondold") == 0)
        tau->psinum = rssringoccs_DR_NewtonDOld;
    else if (strcmp(tau_psitype, "newtondphi") == 0)
        tau->psinum = rssringoccs_DR_NewtonDPhi;
    else if (strcmp(tau_psitype, "newtonperturb") == 0)
        tau->psinum = rssringoccs_DR_NewtonPerturb;
    else if (strcmp(tau_psitype, "ellipse") == 0)
        tau->psinum = rssringoccs_DR_NewtonElliptical;
    else if (strcmp(tau_psitype, "simplefft") == 0)
        tau->psinum = rssringoccs_DR_NewtonSimpleFFT;
    else if (strcmp(tau_psitype, "quadratic") == 0)
        tau->psinum = rssringoccs_DR_NewtonQuadratic;
    else if (strcmp(tau_psitype, "quartic") == 0)
        tau->psinum = rssringoccs_DR_NewtonQuartic;
    else if (strcmp(tau_psitype, "quarticd") == 0)
        tau->psinum = rssringoccs_DR_NewtonDQuartic;
    else if (strcmp(tau_psitype, "fresnel") == 0)
        tau->psinum = rssringoccs_DR_Fresnel;

    /*  strncmp is a C standard library function that compares the first n    *
     *  elements of two strings. If the first seven elements of tau.psitype   *
     *  are "fresnel", but the string is not exactly "fresnel", try to parse  *
     *  the rest of it and extract a value. For example, if                   *
     *  tau.psitype = "fresnel4", try to extract the "4".                     */
    else if (strncmp(tau_psitype, "fresnel", 7) == 0)
    {
        const char *fresnelnum = &tau_psitype[7];
        tau->order = (unsigned char)atol(fresnelnum);
        tau->psinum = rssringoccs_DR_Legendre;

        if (tau->order == 0)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Set_Psitype\n\n"
                "\rCould not parse psitype. atol returned zero. Your input\n"
                "\rhas 'fresnel' in it but either has an invalid entry after,\n"
                "\ror a zero after.\n"
            );
            return;
        }
    }
    else
    {
        char errmes1[1024];
        const char *errmes2 =
            "\r\tquadratic:  Quadratic interpolation of newton-raphson."
            "\r\tcubic:      Cubic interpolation of newton-raphson."
            "\r\tquartic:    Quartic interpolation of newton-raphson."
            "\r\tquarticd:   Quartic interpolation with D perturbation."
            "\r\tellipse:    Newton-Raphson with elliptical perturbation.\n"
            "\r\tfresnel:    Quadratic Fresnel approximation\n"
            "\r\tfresneln:   Legendre polynomial approximation with 1<n<256\n";

        strcpy(
            errmes1,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psitype\n\n"
            "\rIllegal string for psitype. Allowed strings:\n"
            "\r\tnewton:     Newton-Raphson method\n"
            "\r\tnewtond:    Newton-Raphson with D perturbation.\n"
            "\r\tnewtondold: Newton-Raphson with the old D algorithm.\n"
            "\r\tnewtondphi: Newton-Raphson with dD/dphi perturbation.\n"
            "\r\tsimplefft:  A single FFT of the entire data set.\n"
        );

        strcat(errmes1, errmes2);

        tau->order = 0;
        tau->psinum = rssringoccs_DR_None;
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(errmes1);
    }

    free(tau_psitype);
}
