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

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf.h>
#include <rss_ringoccs/include/rss_ringoccs_interpolate.h>

static void write_val(FILE *fp, double x, double y)
{
    fwrite(&x, sizeof(x), 1, fp);
    fwrite(&y, sizeof(y), 1, fp);
}

static double frand(void)
{
    return (double)rand()/RAND_MAX;
}

int main(void)
{
    FILE *linear_fp, *interp_fp;

    unsigned long N = 9;
    unsigned long N_new = 20;
    unsigned long n;

    double start  = -0.5*(N - 1);
    double end    =  0.5*(N - 1);
    double dx     =  (end - start) / N_new;
    double *x     = malloc(sizeof(x) * N);
    double *y     = malloc(sizeof(y) * N);
    double *x_new = malloc(sizeof(x_new) * N_new);
    double *y_new = malloc(sizeof(y_new) * N_new);

    for (n=0; n<N; ++n)
    {
        x[n] = start + n;
        y[n] = frand();
    }

    for (n=0; n<N_new; ++n)
        x_new[n] = start + n*dx;

    rssringoccs_Double_Sorted_Interp1d(x, y, N, x_new, y_new, N_new);

    linear_fp = fopen("linear_binary", "w");
    interp_fp = fopen("interp_binary", "w");

    for (n=0; n<N; ++n)
        write_val(linear_fp, x[n], y[n]);

    for (n=0; n<N_new; ++n)
        write_val(interp_fp, x_new[n], y_new[n]);

    fclose(interp_fp);
    fclose(linear_fp);

    system("graph -T ps -I d < linear_binary > data.ps");
    system("graph -T ps -I d < interp_binary > interp.ps");
    system("rm -f linear_binary interp_binary");

    free(x);
    free(y);
    free(x_new);
    free(y_new);

    return 0;
}
