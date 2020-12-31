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
 ******************************************************************************/

#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <stdio.h>
#include <stdlib.h>

static void write_val(FILE *fp, double x, double y)
{
    fwrite(&x, sizeof(x), 1, fp);
    fwrite(&y, sizeof(y), 1, fp);
}

int main(void)
{
    rssringoccs_DLPCSV *dlp;
    unsigned long n, start_n, end_n;
    double start = 87400.0;
    double end   = 87600.0;
    FILE *fp;

    dlp = rssringoccs_Get_DLP("../Test_Data/Rev007E_X43_Maxwell_DLP_500M.TAB",
                              rssringoccs_False);
    if (dlp == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\ttest_get_dlp_csv\n\n"
             "rssringoccs_Get_DLP returned NULL.\n");
        return -1;
    }
    else if (dlp->error_occurred)
    {
        if (dlp->error_message == NULL)
        {
            puts("Error Encountered: rss_ringoccs\n"
                "\ttest_get_dlp_csv\n\n"
                "dlp->error_occurred set to true with no error message.\n");
            rssringoccs_Destroy_DLPCSV(&dlp);
            return -1;
        }
        else
        {
            printf("Error Encountered: rss_ringoccs\n"
                   "\ttest_get_dlp_csv\n\n"
                   "dlp->error_occurred set to true. Printing error:\n\n%s",
                   dlp->error_message);
            rssringoccs_Destroy_DLPCSV(&dlp);
            return -1;
        }
    }

    fp = fopen("rev007_dlp_binary", "w");

    n = 0;
    while (dlp->rho_km_vals[n] < start)
        n++;

    start_n = n;

    while(dlp->rho_km_vals[n] < end)
        n++;

    end_n = n;

    for (n=start_n; n<end_n; ++n)
        write_val(fp, dlp->rho_km_vals[n], dlp->p_norm_vals[n]);

    fclose(fp);

    system("graph --font-size 0.03 -T ps -I d < rev007_dlp_binary > plot.ps");
    system("rm -f rev007_dlp_binary");
    rssringoccs_Destroy_DLPCSV(&dlp);
    return 0;
}
