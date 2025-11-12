#include <libtmpl/include/types/tmpl_cyl_fresnel_geometry_double.h>
#include <libtmpl/include/types/tmpl_complex_double.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

tmpl_ComplexDouble
rssringoccs_Fresnel_Kernel(const rssringoccs_TAUObj * const tau,
                           size_t offset,
                           size_t center)
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
