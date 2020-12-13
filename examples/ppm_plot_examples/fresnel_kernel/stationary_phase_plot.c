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

#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>
#include <stdio.h>

#define CYAN (unsigned char)204, (unsigned char)238, (unsigned char)255
#define GREEN (unsigned char)204, (unsigned char)255, (unsigned char)153
#define BLACK (unsigned char)0, (unsigned char)0, (unsigned char)0

int main(void)
{
    FILE *fp0, *fp1;

    /*  Geometry values from Rev007 E.                                        */
    double r0   = 88734.619072;
    double phi0 = 251.115963 * rssringoccs_One_Pi / 180.0;
    double B    = -23.571288 * rssringoccs_One_Pi / 180.0;
    double x0   = r0 * rssringoccs_Double_Cos(phi0);
    double y0   = r0 * rssringoccs_Double_Sin(phi0);
    double xc   = -291898.181305;
    double yc   = -83967.945080;
    double zc   =  114820.549365;
    double diff = 100.0;

    double x_min = x0 - diff;
    double x_max = x0 + diff;
    double y_min = y0 - diff;
    double y_max = y0 + diff;
    double xs, ys;
    double dpsi, dpsi0;
    double r, phi;
    double err, err0;
    double EPS = 10.0;

    unsigned int size = 1024;
    unsigned int x, y;

    double ds         = diff/size;
    double rcp_factor = 1.0/(size-1.0);

    double dist;
    double D = rssringoccs_Double_Sqrt((xc-x0)*(xc-x0) +
                                       (yc-y0)*(yc-y0) +
                                       zc*zc);

    double freq = 8426334772.410973;
    double lambda = rssringoccs_Speed_Of_Light_KMS / freq;

    double kD;
    double kD0 = rssringoccs_Two_Pi * D / lambda;

    fp0 = fopen("stationary_phase_plot0.ppm", "w");
    fprintf(fp0, "P6\n%d %d\n255\n", size, size);

    fp1 = fopen("stationary_phase_plot1.ppm", "w");
    fprintf(fp1, "P6\n%d %d\n255\n", size, size);

    /*  Loop through each pixel.                                              */
    for (y=0; y<size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        ys = (y - 1.0) * (y_max - y_min) * rcp_factor + y_min;

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            xs = x * (x_max - x_min) * rcp_factor + x_min;

            phi = rssringoccs_Double_Arctan2(ys, xs);
            r   = rssringoccs_Double_Sqrt(xs*xs + ys*ys);

            /*  Compute the complex number z_x + i z_y.                       */
            dist = rssringoccs_Double_Sqrt((xc-xs)*(xc-xs) +
                                           (yc-ys)*(yc-ys) +
                                           zc*zc);

            kD   = rssringoccs_Two_Pi * dist / lambda;
            dpsi = rssringoccs_Double_Fresnel_dPsi_dPhi(kD, r, r0, phi,
                                                         phi0, B, dist);

            dpsi0 = rssringoccs_Double_Fresnel_dPsi_dPhi(kD0, r, r0, phi,
                                                         phi0, B, D);
            err  = rssringoccs_Double_Abs(dpsi);
            err0 = rssringoccs_Double_Abs(dpsi0);

            rssringoccs_RGB_Linear_Gradient(err0, 0.0, 5.0e5, fp0);
            rssringoccs_RGB_Linear_Gradient(err,  0.0, 5.0e5, fp1);
        }
        printf("%u %f %f\n", y, err, err0);
    }
    return 0;
}
