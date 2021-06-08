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

double rssringoccs_Double_Fresnel_Psi_Old(double kD, double r, double r0,
                                          double phi, double phi0,
                                          double B, double D)
{
    double xi, eta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrt(1.0+eta-2.0*xi) + xi - 1.0);
}
