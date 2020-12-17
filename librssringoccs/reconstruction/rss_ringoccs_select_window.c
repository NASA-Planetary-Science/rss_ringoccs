#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      select_window_func                                                    *
 *  Purpose:                                                                  *
 *      Parse the wtype variable of a DLPObj pointer and retrieve the         *
 *      corresponding function.                                               *
 *  Arguments:                                                                *
 *      fw (rss_ringoccs_window_func):                                        *
 *          Function pointer for the window function.                         *
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
void select_window_func(rss_ringoccs_window_func *fw,
                        rssringoccs_TAUObj *tau)
{
    /*  First, check that wtype is even valid! If dlp->wtype is a NULL        *
     *  pointer and we try to strcmp it, we'll get a segmentation fault.      *
     *  Check that it's not null with "if (!dlp->wtype)"                      */
    if (!tau->wtype)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tselect_window_func\n\n"
            "\r\ttau->wtype is NULL. Returning.\n\n",
            tau
        );

    /*  If dlp->wtype is not null, check if it's a valid string and cast the  *
     *  selected window type to the fw pointer.                               */
    else if (strcmp(tau->wtype, "rect") == 0)
        *fw = rssringoccs_Double_Rect_Window;
    else if (strcmp(tau->wtype, "coss") == 0)
        *fw = rssringoccs_Double_Coss_Window;
    else if (strcmp(tau->wtype, "kb20") == 0)
        *fw = rssringoccs_Double_Kaiser_Bessel_2_0;
    else if (strcmp(tau->wtype, "kb25") == 0)
        *fw = rssringoccs_Double_Kaiser_Bessel_2_5;
    else if (strcmp(tau->wtype, "kb35") == 0)
        *fw = rssringoccs_Double_Kaiser_Bessel_3_5;
    else if (strcmp(tau->wtype, "kbmd20") == 0)
        *fw = rssringoccs_Double_Modified_Kaiser_Bessel_2_0;
    else if (strcmp(tau->wtype, "kbmd25") == 0)
        *fw = rssringoccs_Double_Modified_Kaiser_Bessel_2_5;
    else if (strcmp(tau->wtype, "kbmd35") == 0)
        *fw = rssringoccs_Double_Modified_Kaiser_Bessel_3_5;
    else
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tselect_window_func\n\n"
            "\r\tInvalid input for tau->wtype. Returning.\n\n",
            tau
        );
}
