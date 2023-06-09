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
#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdio.h>
#include <stdlib.h>

static void write_val(FILE *fp, double x, double y)
{
    fwrite(&x, sizeof(x), 1, fp);
    fwrite(&y, sizeof(y), 1, fp);
}

int main(void)
{
    const char *geo = "../Test_Data/Rev007E_X43_Maxwell_GEO.TAB";
    const char *cal = "../Test_Data/Rev007E_X43_Maxwell_CAL.TAB";
    const char *dlp = "../Test_Data/Rev007E_X43_Maxwell_DLP_500M.TAB";

    rssringoccs_CSVData *data;
    size_t n;
    FILE *fp;

    data = rssringoccs_Extract_CSV_Data(geo, cal, dlp, NULL, tmpl_False);

    if (data == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\ttest_extract_csv_data\n\n"
             "rssringoccs_Extract_CSV_Data returned NULL.\n");
        return -1;
    }

    if (data->error_occurred)
    {
        if (data->error_message == NULL)
        {
            puts("Error Encountered: rss_ringoccs\n"
                "\ttest_extract_csv_data\n\n"
                "data->error_occurred set to true with no error message.\n");
            rssringoccs_Destroy_CSV(&data);
            return -1;
        }

        printf("Error Encountered: rss_ringoccs\n"
               "\ttest_extract_csv_data\n\n"
               "data->error_occurred set to true. Printing error:\n\n%s",
               data->error_message);

        rssringoccs_Destroy_CSV(&data);
        return -1;
    }

    fp = fopen("rev007_plot_binary", "w");

    for (n = (size_t)0; n < data->n_elements; ++n)
        write_val(fp, data->rho_km_vals[n], data->D_km_vals[n]);

    fclose(fp);

    system("graph --font-size 0.03 -T ps -I d < rev007_plot_binary > plot.ps");
    system("rm -f rev007_plot_binary");
    rssringoccs_Destroy_CSV(&data);
    return 0;
}
