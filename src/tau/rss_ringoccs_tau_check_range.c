/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Is_NaN and tmpl_Double_Is_Inf declared here.                  */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for checking if the input range is valid.                        */
void rssringoccs_Tau_Check_Range(rssringoccs_TAUObj * const tau)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!tau)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (tau->error_occurred)
        return;

    /*  Range values should be real numbers. Check for NaN (Not-a-Number).    */
    if (tmpl_Double_Is_NaN(tau->range[0]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[0] is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Same check for the next value in the range list.                      */
    if (tmpl_Double_Is_NaN(tau->range[1]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[1] is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  The range should also be finite. Treat infinity as an error.          */
    if (tmpl_Double_Is_Inf(tau->range[0]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[0] is infinite.\n\n";

        return;
    }

    /*  Same check for the next value in the range list.                      */
    if (tmpl_Double_Is_Inf(tau->range[1]))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rrange[1] is infinite.\n\n";

        return;
    }

    /*  The lower range should be positive.                                   */
    if (tau->range[0] < 0.0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rStarting value for range is negative.\n\n";

        return;
    }

    /*  Lastly, the range list should be in increasing order.                 */
    if (tau->range[0] > tau->range[1])
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Range\n\n"
            "\rStarting value for range is greater than final value.\n\n";

        return;
    }
}
/*  End of rssringoccs_Tau_Check_Range.                                       */
