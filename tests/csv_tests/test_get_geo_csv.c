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
    rssringoccs_GeoCSV *geo;
    unsigned long n;
    FILE *fp;

    geo = rssringoccs_Get_Geo(
        "../Test_Data/Rev007E_X43_Maxwell_GEO.TAB", tmpl_False
    );

    if (geo == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\ttest_get_geo_csv\n\n"
             "rssringoccs_Get_Geo returned NULL.\n");
        return -1;
    }
    else if (geo->error_occurred)
    {
        if (geo->error_message == NULL)
        {
            puts("Error Encountered: rss_ringoccs\n"
                "\ttest_get_geo_csv\n\n"
                "geo->error_occurred set to true with no error message.\n");
            rssringoccs_Destroy_GeoCSV(&geo);
            return -1;
        }
        else
        {
            printf("Error Encountered: rss_ringoccs\n"
                   "\ttest_get_geo_csv\n\n"
                   "geo->error_occurred set to true. Printing error:\n\n%s",
                   geo->error_message);
            rssringoccs_Destroy_GeoCSV(&geo);
            return -1;
        }
    }

    fp = fopen("rev007_plot_binary", "w");

    for (n=0; n<geo->n_elements; ++n)
        write_val(fp, geo->rho_km_vals[n], geo->F_km_vals[n]);

    fclose(fp);

    system("graph --font-size 0.03 -T ps -I d < rev007_plot_binary > plot.ps");
    system("rm -f rev007_plot_binary");
    rssringoccs_Destroy_GeoCSV(&geo);
    return 0;

}
