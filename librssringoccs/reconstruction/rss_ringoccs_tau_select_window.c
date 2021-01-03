#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      select_window_func                                                    *
 *  Purpose:                                                                  *
 *      Parse the wtype variable of a DLPObj pointer and retrieve the         *
 *      corresponding function.                                               *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          rss_ringoccs_diffraction_correction.h.                            *
 *  Notes:                                                                    *
 *      1.) This is a "static" function meanings its use is limited to this   *
 *          file. Its only real purpose is for checking DLP instances which   *
 *          are the inputs to the various diffraction correction functions.   *
 *      2.) If the function fails, dlp->status is set to 5. This means an     *
 *          invalid string for wtype was given.
 ******************************************************************************/
void rssringoccs_Tau_Select_Window_Func(rssringoccs_TAUObj *tau)
{
    if (!tau)
        return;

    /*  First, check that wtype is even valid! If dlp->wtype is a NULL        *
     *  pointer and we try to strcmp it, we'll get a segmentation fault.      *
     *  Check that it's not null with "if (!dlp->wtype)"                      */
    if (!tau->wtype)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tselect_window_func\n\n"
            "\r\ttau->wtype is NULL. Returning.\n\n"
        );
    }

    /*  If dlp->wtype is not null, check if it's a valid string and cast the  *
     *  selected window type to the window_func pointer.                      */
    else if (strcmp(tau->wtype, "rect") == 0)
        tau->window_func = rssringoccs_Double_Rect_Window;
    else if (strcmp(tau->wtype, "coss") == 0)
        tau->window_func = rssringoccs_Double_Coss_Window;
    else if (strcmp(tau->wtype, "kb20") == 0)
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_2_0;
    else if (strcmp(tau->wtype, "kb25") == 0)
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_2_5;
    else if (strcmp(tau->wtype, "kb35") == 0)
        tau->window_func = rssringoccs_Double_Kaiser_Bessel_3_5;
    else if (strcmp(tau->wtype, "kbmd20") == 0)
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_2_0;
    else if (strcmp(tau->wtype, "kbmd25") == 0)
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_2_5;
    else if (strcmp(tau->wtype, "kbmd35") == 0)
        tau->window_func = rssringoccs_Double_Modified_Kaiser_Bessel_3_5;
    else
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tselect_window_func\n\n"
            "\r\tInvalid input for tau->wtype. Returning.\n\n"
        );
    }
}
