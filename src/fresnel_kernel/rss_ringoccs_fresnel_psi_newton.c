/******************************************************************************
 *                                 LICENSE                                    *
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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 8, 2021                                                  *
 ******************************************************************************/

#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double
rssringoccs_Newton_Raphson_Fresnel_Psi(double k, double r, double r0,
                                       double phi, double phi0, double B,
                                       double D, double EPS, unsigned int toler)
{
    double dphi, err;
    unsigned int n = 0U;
    dphi  = rssringoccs_Fresnel_dPsi_dPhi(k, r, r0, phi, phi0, B, D);
    err = fabs(dphi);
    while(err > EPS)
    {
        dphi = rssringoccs_Fresnel_dPsi_dPhi(k, r, r0, phi, phi0, B, D);
        phi -= dphi/rssringoccs_Fresnel_d2Psi_dPhi2(k, r, r0, phi, phi0, B, D);
        ++n;
        if (n > toler)
            break;

        err = fabs(dphi);
    }
    return phi;
}
