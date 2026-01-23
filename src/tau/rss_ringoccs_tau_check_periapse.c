/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Is_NaN and tmpl_Double_Is_Inf declared here.                  */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking if the input periapse is valid.                     */
void rssringoccs_Tau_Check_Periapse(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  The periapse should be a real number. Check for NaN (Not-a-Number).   */
    if (tmpl_Double_Is_NaN(tau->periapse))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Periapse\n\n"
            "\rperiapse is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Periapse should also be finite. Treat infinity as an error.           */
    if (tmpl_Double_Is_Inf(tau->periapse))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Periapse\n\n"
            "\rperiapse is infinite.\n\n";

        return;
    }

    /*  The periapse is allowed to be between -2pi and 2pi, inclusive.        */
    if (tau->periapse < -tmpl_double_two_pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse less than -2pi.\n\n";

        return;
    }

    /*  Same check, peripase should be bounded by +2 pi.                      */
    if (tau->periapse > tmpl_double_two_pi)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Keywords\n\n"
            "\rInput periapse greater than 2pi.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Periapse.                                    */
