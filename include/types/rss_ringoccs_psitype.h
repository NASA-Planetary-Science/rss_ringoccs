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
 *  Purpose:                                                                  *
 *      Provides an enum listing all of the types of Fresnel inversion.       *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       December 28, 2025                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_PSITYPE_H
#define RSS_RINGOCCS_TYPES_PSITYPE_H

/*  All of the types of reconstruction.                                       */
enum rssringoccs_PsiType {

    /*  Fresnel quadratic approximation.                                      */
    rssringoccs_PsiType_Fresnel,

    /*  Legendre polynomials, higher order Fresnel approximations.            */
    rssringoccs_PsiType_Legendre,

    /*  Newton-Raphson, using a single FFT across the entire data set.        */
    rssringoccs_PsiType_NewtonSimpleFFT,

    /*  Newton-Raphson method, slow but accurate, no interpolation performed. */
    rssringoccs_PsiType_NewtonRiemann,

    /*  Newton-Raphson with a linear Filon-like quadrature method.            */
    rssringoccs_PsiType_NewtonLinearFilon,

    /*  Newton-Raphson with a quadratic Filon-like quadrature method.         */
    rssringoccs_PsiType_NewtonQuadraticFilon,

    /*  Newton-Raphson with elliptical corrections, and interpolations.       */
    rssringoccs_PsiType_EllipticNewton,

    /*  Newton-Raphson with arbitrary quartic perturbation polynomial.        */
    rssringoccs_PsiType_NewtonPerturb,

    /*  Even degree interpolations for the Newton-Raphson method.             */
    rssringoccs_PsiType_Newton4,
    rssringoccs_PsiType_Newton8,
    rssringoccs_PsiType_Newton16,

    /*  Even degree interpolations for the elliptic Newton-Raphson method.    */
    rssringoccs_PsiType_EllipticNewton4,
    rssringoccs_PsiType_EllipticNewton8,
    rssringoccs_PsiType_EllipticNewton16,

    /*  Indicates an error.                                                   */
    rssringoccs_PsiType_None
};

#endif
/*  End of include guard.                                                     */
