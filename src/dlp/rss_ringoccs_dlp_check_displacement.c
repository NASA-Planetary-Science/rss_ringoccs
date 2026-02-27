/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *                     rss_ringoccs_dlp_check_displacement                    *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks for errors in the displacement found in a dlp object.          *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_DLP_Check_Displacement                                    *
 *  Purpose:                                                                  *
 *      Checks for a few common errors found in the displacement.             *
 *  Arguments:                                                                *
 *      dlp (rssringoccs_DLPObj * const):                                     *
 *          A pointer to a dlp object.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_math.h:                                                          *
 *          tmpl_Double_Abs:                                                  *
 *              Computes the absolute value of a real number.                 *
 *          tmpl_Double_Is_Inf:                                               *
 *              Checks if a double is +/- infinity.                           *
 *          tmpl_Double_Is_NaN:                                               *
 *              Checks if a double is NaN (Not-a-Number).                     *
 *  Notes:                                                                    *
 *      1.) If the error_occurred Boolean was previously set to true, this    *
 *          function does nothing and skips all checks.                       *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          A detailed description of the geometry of a ring occultation can  *
 *          be found here, and this includes the definition of rho.           *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_math.h:                                                          *
 *          tmpl_Double_Abs declared here.                                    *
 *  3.) rss_ringoccs_dlp.h:                                                   *
 *          DLP definition and prototype for the function given here.         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       April 11, 2025                                                *
 ******************************************************************************/

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_Double_Abs declared here.                                            */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Checks the displacement in a dlp object for errors.                       */
void rssringoccs_DLP_Check_Displacement(rssringoccs_DLPObj * const dlp)
{
    /*  If the input is NULL there is nothing to be done.                     */
    if (!dlp)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (dlp->error_occurred)
        return;

    /*  The displacement should be real. Check for NaN (Not-a-Number).        */
    if (tmpl_Double_Is_NaN(dlp->dx_km))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Displacement\n\n"
            "\rdx_km is NaN (Not-a-Number).\n\n";

        return;
    }

    /*  Displacement should also be finite. Treat infinity as an error.       */
    if (tmpl_Double_Is_Inf(dlp->dx_km))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Displacement\n\n"
            "\rdx_km is infinite.\n\n";

        return;
    }

    /*  The sample spacing needs to be non-zero. For one, there is a finite   *
     *  amount of data available, and secondly we divide by dx_km to compute  *
     *  the length of the window array. Avoid divide-by-zero, check for this. */
    if (dlp->dx_km == 0.0)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_DLP_Check_Displacement\n\n"
            "\rdx_km = rho_km_vals[1] - rho_km_vals[0] = 0.\n\n";

        return;
    }
}
/*  End of rssringoccs_DLP_Check_Displacement.                                */
