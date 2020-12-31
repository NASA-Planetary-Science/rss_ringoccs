/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *  Purpose:                                                                  *
 *      Provide tests for the accuracy and efficiency of rss_ringoccs         *
 *      absolute value function compared to the one provided by the C99       *
 *      standard in math.h.                                                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <gsl/gsl_sf.h>

/*  The comparison functions are found here.                                  *
 *  NOTE:                                                                     *
 *      You will need to compile librssringoccs_compare to run this. This can *
 *      be done via the config_librssringoccs_compare.sh script found in      *
 *      tests/librssringoccs_compare. We can then link this library via       *
 *      -lrssringoccs_compare (see below).                                    */
#include <rss_ringoccs/include/rss_ringoccs_compare_funcs.h>

int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    double start = -1000;
    double end   =  1000;

    /*  We'll test on 100 million points between start and end.               */
    unsigned long N = 1e8;

    /*  Use the compare function to test rssringoccs_Double_Abs against fabs. */
    rssringoccs_Compare_Double_Funcs("gsl-2.6",
                                     gsl_sf_bessel_J0,
                                     "rss_ringoccs",
                                     rssringoccs_Double_Bessel_J0,
                                     start, end, N);

    return 0;
}
/*  End of main.                                                              */
