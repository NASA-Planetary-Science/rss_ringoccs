#include <stdlib.h>

/*  To help manipulate strings, both native C strings const char* and char*,  *
 *  as well as Python strings passed by the user, include string.h. This is a *
 *  standard library header containing strcmp (string-compare), strncmp       *
 *  (compares first n characters of string), strcpy (string-copy), strncpy,   *
 *  and strlen (string-length), all of which are used at some point.          */
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void
rssringoccs_Tau_Set_Psitype(const char *psitype, rssringoccs_TAUObj* tau)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (tau->psitype != NULL)
        free(tau->psitype);

    tau->psitype = rssringoccs_strdup(psitype);
    rssringoccs_Remove_Spaces(tau->psitype);
    rssringoccs_Make_Lower(tau->psitype);

    if (strcmp(tau->psitype, "newton") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_Newton;
    }
    else if (strcmp(tau->psitype, "newtond") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_NewtonD;
    }
    else if (strcmp(tau->psitype, "newtondold") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_NewtonDOld;
    }
    else if (strcmp(tau->psitype, "newtondphi") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_NewtonDPhi;
    }
    else if (strcmp(tau->psitype, "ellipse") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_Elliptical;
    }

    /*  If psitype is just "fresnel", set order to 1. The code will           *
     *  automatically pass the dlp instance to the correct function.          */
    else if (strcmp(tau->psitype, "fresnel") == 0)
    {
        tau->order  = 0;
        tau->psinum = rssringoccs_DR_Fresnel;
    }

    /*  strncmp is a C standard library function that compares the first n    *
     *  elements of two strings. If the first seven elements of tau.psitype   *
     *  are "fresnel", but the string is not exactly "fresnel", try to parse  *
     *  the rest of it and extract a value. For example, if                   *
     *  tau.psitype = "fresnel4", try to extract the "4".                     */
    else if (strncmp(tau->psitype, "fresnel", 7) == 0)
    {
        const char *fresnelnum = &tau->psitype[7];
        tau->order = strtol(fresnelnum, NULL, 10);
        tau->psinum = rssringoccs_DR_Legendre;
    }
    else
    {
        tau->order = 0;
        tau->psinum = rssringoccs_DR_None;
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Set_Psitype\n\n"
            "\r\tIllegal string for psitype. Allowed strings are:\n"
            "\r\t\tnewton:     Newton-Raphson method\n"
            "\r\t\tnewtond:    Newton-Raphson with D perturbation.\n"
            "\r\t\tnewtondold: Newton-Raphson with the old D algorithm.\n"
            "\r\t\tnewtondphi: Newton-Raphson with dD/dphi perturbation.\n"
            "\r\t\tellipse:    Newton-Raphson with elliptical perturbation.\n"
            "\r\t\tfresnel:    Quadratic Fresnel approximation\n"
            "\r\t\tfresneln:   Legendre polynomial approximation with 1<n<256\n"
        );
    }
}
