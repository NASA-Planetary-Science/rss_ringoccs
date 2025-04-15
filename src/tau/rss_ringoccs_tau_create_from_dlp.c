/*  NULL pointers are given here.                                             */
#include <stdlib.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for creating a Tau object from DLP data and a resolution.        */
rssringoccs_TAUObj *
rssringoccs_Tau_Create_From_DLP(const rssringoccs_DLPObj *dlp, double res)
{
    /*  Try to allocate memory for a new tau object.                          */
    rssringoccs_TAUObj *tau = malloc(sizeof(*tau));

    /*  malloc returns NULL on failure. Check for this.                       */
    if (!tau)
        return NULL;

    /*  Initialize all members to their zero values or their defaults.        */
    rssringoccs_Tau_Init(tau);

    /*  Check if the input dlp is NULL. This is treated as an error.          */
    if (!dlp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Create_From_DLP\n\n"
            "\rInput dlp is NULL.\n\n";

        return tau;
    }

    /*  Check if the input dlp has its error_occurred member set to true.     *
     *  Again, this is treated as an error for the Tau object.                */
    if (dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Create_From_DLP\n\n"
            "\rInput dlp has error_occurred set to true.\n\n";

        return tau;
    }

    /*  Check that the resolution is a legal value.                           */
    if (res <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Create_From_DLP\n\n"
            "\rInput resolution is not positive.\n\n";

        return tau;
    }

    /*  Set the resolution variable to the user-provided input.               */
    tau->resolution_km = res;

    /*  Most of the variables for the Tau object can be computed directly     *
     *  from the data in the DLP object. The only variables that are not set  *
     *  after this function is called (and are still their zero values or     *
     *  NULL) are the reconstruction and forward modeling variables.          */
    rssringoccs_Tau_Copy_DLP_Data(dlp, tau);

    /*  Check the data for possible errors before returning.                  */
    rssringoccs_Tau_Check_Geometry(tau);
    return tau;
}
/*  End of rssringoccs_Tau_Create_From_DLP.                                   */
