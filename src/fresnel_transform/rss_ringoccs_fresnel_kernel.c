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
 *                         rss_ringoccs_fresnel_kernel                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Computes the stationary Fresnel kernel for a given point in the data. *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Fresnel_Kernel                                            *
 *  Purpose:                                                                  *
 *      Computes the stationary Fresnel kernel from a given Tau object.       *
 *  Arguments:                                                                *
 *      tau (const rssringoccs_TAUObj * const tau):                           *
 *          The Tau object with all of the geometry and diffraction data.     *
 *      center (const size_t):                                                *
 *          The index for the center of the window currently being considered.*
 *          For the forward transform this is the index for "rho0", and for   *
 *          the inverse transform this is the index for "rho".                *
 *      offset (const size_t):                                                *
 *          The index for the variable of integration which runs across the   *
 *          window. This is "rho" for the forward transform and "rho0" for    *
 *          the inverse transform.                                            *
 *  Output:                                                                   *
 *      kernel (tmpl_ComplexDouble):                                          *
 *          The stationary Fresnel kernel.                                    *
 *  Called Functions:                                                         *
 *      tmpl_complex.h:                                                       *
 *          tmpl_CDouble_ConjugateSelf:                                       *
 *              Conjugates a complex number, storing the result in the input. *
 *      tmpl_cyl_fresnel_optics.h:                                            *
 *          tmpl_Double_Stationary_Cyl_Fresnel_Kernel:                        *
 *              Computes the stationary Fresnel kernel using Newton's method. *
 *      tmpl_vec2.h:                                                          *
 *          tmpl_2DDouble_Polard:                                             *
 *              Creates a vector from polar coordinates with angle in degrees.*
 *      tmpl_vec3.h:                                                          *
 *          tmpl_3DDouble_Rect:                                               *
 *              Returns the vector (x, y, z) given three real numbers.        *
 *  Method:                                                                   *
 *      Using the geometry data found in the Tau object, construct a          *
 *      tmpl_CylFresnelGeometryDouble object that can then be passed to       *
 *      libtmpl for the main computation. If the forward transform is being   *
 *      used, then conjugate the end result.                                  *
 *  Notes:                                                                    *
 *      1.) There are no checks for NULL points. It is assumed tau is valid.  *
 *                                                                            *
 *      2.) There are no checks for errors. This function is called inside    *
 *          the main for-loop of the various transforms. It is the users      *
 *          obligation to ensure the input data is valid.                     *
 *                                                                            *
 *      3.) There are no checks that the center and offset arguments          *
 *          correspond to indices for the actual data.                        *
 *                                                                            *
 *      4.) It is assumed that the input Tau object has all of its angles in  *
 *          degrees. This should always be the case, regardless.              *
 *  References:                                                               *
 *      1.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. The Fresnel kernel is described here.     *
 *                                                                            *
 *      2.) Goodman, J. (2005)                                                *
 *          Introduction to Fourier Optics                                    *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *                                                                            *
 *          Covers most of the theory behind diffraction and the application  *
 *          of Fourier analysis to optics. The Fresnel transform is given an  *
 *          in-depth treatise in this book.                                   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_cyl_fresnel_geometry_double.h:                                   *
 *          Header file providing the tmpl_CylFresnelGeometryDouble typedef.  *
 *  2.) tmpl_cyl_fresnel_optics.h:                                            *
 *          Header file providing the stationary Fresnel kernel function.     *
 *  3.) tmpl_complex.h:                                                       *
 *          Header providing complex numbers types and functions.             *
 *  4.) tmpl_vec2.h:                                                          *
 *          Header file with 2D vectors and vector functions.                 *
 *  4.) tmpl_vec2.h:                                                          *
 *          Header file providing 3D vector types and routines.               *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       November 12, 2025                                             *
 ******************************************************************************/

/*  Geometry object with rho, rho, and the spacecraft position vector.        */
#include <libtmpl/include/types/tmpl_cyl_fresnel_geometry_double.h>

/*  Function for computing the stationary Fresnel phase provided here.        */
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>

/*  Complex number types and functions provided here.                         */
#include <libtmpl/include/tmpl_complex.h>

/*  2D vector types and functions found here.                                 */
#include <libtmpl/include/tmpl_vec2.h>

/*  3D vector types and functions given here.                                 */
#include <libtmpl/include/tmpl_vec3.h>


#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

tmpl_ComplexDouble
rssringoccs_Fresnel_Kernel(const rssringoccs_TAUObj * const tau,
                           const size_t center,
                           const size_t offset)
{
    tmpl_ComplexDouble ker;
    tmpl_CylFresnelGeometryDouble geo;
    const size_t val0 = (tau->use_fwd ? center : offset);
    const size_t val1 = (tau->use_fwd ? offset : center);

    /*  This function does not assume the ideal geometry present in MTR86     *
     *  where u . y = 0, u being the vector from the spacecraft to the ring   *
     *  intercept point. Instead we compute using the full Fresnel kernel.    *
     *  This requires each vector in their Cartesian coordinates.             */
    geo.position = tmpl_3DDouble_Rect(
        tau->rx_km_vals[val0],
        tau->ry_km_vals[val0],
        tau->rz_km_vals[val0]
    );

    geo.intercept = tmpl_2DDouble_Polard(
        tau->rho_km_vals[val0], tau->phi_deg_vals[val0]
    );

    geo.dummy = tmpl_2DDouble_Polard(
        tau->rho_km_vals[val1], tau->phi_deg_vals[val0]
    );

    /*  The full Fresnel kernel is given by the stationary phase method:      *
     *                                                                        *
     *                                                                        *
     *        2 pi                                                            *
     *         -                          _________                           *
     *        | | exp(i  psi)            / 2 pi     exp(i (psi_s - pi/4))     *
     *        |   ----------- dphi ~=   / --------- -------------------       *
     *      | |   | R - rho |         \/  |psi_s''|    | R - rho_s |          *
     *       -                                                                *
     *       0                                                                *
     *                                                                        *
     *  libtmpl has tools to compute this using Newton's Method. Return this. */
    ker = tmpl_Double_Stationary_Cyl_Fresnel_Kernel(
        tau->k_vals[val0],
        &geo,
        tau->EPS,
        tau->toler
    );

    /*  The inverse transform is approximated by negating the Fresnel phase:  *
     *                                                                        *
     *                   -             _________                              *
     *                  | | ^         / 2 pi     exp(-i (psi_s-pi/4))         *
     *      T(rho) ~=   |   T(r0) r  / --------- -------------------- dr0     *
     *                | |          \/  |psi_s''|     | R - rho_s |            *
     *                 -                                                      *
     *                                                                        *
     *  The kernel for the inverse transform is hence the complex conjugate   *
     *  of the kernel for the forward transform.                              */
    if (!tau->use_fwd)
        tmpl_CDouble_ConjugateSelf(&ker);

    return ker;
}
