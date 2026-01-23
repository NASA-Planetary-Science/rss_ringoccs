/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Is_NaN and tmpl_Double_Is_Inf declared here.                  */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking if the input Allan deviation is valid.              */
void rssringoccs_Tau_Check_Allan_Deviation(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  Sigma should be a real number. Check for NaN (Not-a-Number).          */
    if (tmpl_Double_Is_NaN(tau->sigma))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Allan_Deviation\n\n"
            "\rAllan deviation is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Sigma should also be finite. Treat infinity as an error.              */
    if (tmpl_Double_Is_Inf(tau->sigma))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Allan_Deviation\n\n"
            "\rAllan deviation is infinite.\n\n";

        return;
    }

    /*  The Allen deviation must be positive.                                 */
    if (tau->sigma <= 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Allan_Deviation\n\n"
            "\rInput sigma (Allen deviation) is not positive.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Allan_Deviation.                             */
