#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      check_data_range                                                      *
 *  Purpose:                                                                  *
 *      Check if the window ranges needed for reconstruction are permitted by *
 *      the available data. This check is important to avoid segmentation     *
 *      faults at runtime.                                                    *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          rss_ringoccs_diffraction_correction.h.                            *
 *  Output:                                                                   *
 *      Out (Boolean):                                                        *
 *          An integer, either 0 or 1 depending on if the checks passed.      *
 *  Notes:                                                                    *
 *      1.) This is a "static" function meanings its use is limited to this   *
 *          file. Its only real purpose is for checking DLP instances which   *
 *          are the inputs to the various diffraction correction functions.   *
 *      2.) If the check fails, there is something wrong with the input data. *
 *          Removing this check may result in a "segmentation fault 11" error *
 *          at runtime with no additional error message, a scenario one would *
 *          like to avoid.                                                    *
 ******************************************************************************/
void check_tau_data_range(rssringoccs_TAUObj *tau, double two_dx)
{
    /* Create variables to check what the minimum and maximum indices are.    */
    long current_min, current_max;

    /* Variable for keeping track of the largest window needed in the loop.   */
    long max_nw_pts = 0;

    /* Variables to keep track of the min and max window indices.             */
    long min_requested = tau->start;
    long max_requested = min_requested + tau->n_used;

    /* Variable for indexing the for loop.                                    */
    long n;

    /* Variables for the window size and number of points in the window.      */
    double win_size;
    long nw_pts;

    /* Loop through every point, check window width, and ensure you have      *
     * enough data to the left and right for data processing.                 */
    for (n=tau->start; n<=max_requested; ++n)
    {
        /*  Gather the window size and number of points in the window.        */
        win_size = tau->w_km_vals[n];
        nw_pts   = ((long)(win_size / two_dx))+1;

        /*  The smallest radius needed for the window corresponds to the      *
         *  current point minus the number of points in the window. Similarly *
         *  the maximum corresponds to the sum. Compute these.                */
        current_min = n-nw_pts;
        current_max = n+nw_pts;

        /*  If current_min is smaller than min_requested reset min_requested  *
         *  and similarly if current_max is larger than max_requested.        */
        if (current_min < min_requested)
            min_requested   = current_min;

        if (current_max < max_requested)
            max_requested   = current_max;

        /*  If nw_pts got larger, record this in max_nw_pts.                  */
        if (nw_pts > max_nw_pts)
            max_nw_pts = nw_pts;
    }

    /*  If min_requested is negative, the window is too large. Similarly, if  *
     *  max_requested goes between the size of the array.                     */
    if (min_requested < 0)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data_range\n\n"
            "\rRequested data range goes beyond the range available. The\n"
            "\rrequested minimum in radius is less than the minimum radius\n"
            "\ravailable. Returning.\n\n",
            tau
        );
    else if (max_requested > (long)tau->arr_size)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data_range\n\n"
            "\rRequested data range goes beyond the range available. The\n"
            "\rrequested maximum in radius is greater than the maximum radius\n"
            "\ravailable. Returning.\n\n",
            tau
        );
}
