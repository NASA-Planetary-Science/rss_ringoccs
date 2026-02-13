#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/types/tmpl_cyl_fresnel_geometry_double.h>
#include <libtmpl/include/types/tmpl_complex_double.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

double
rssringoccs_Fast_Fresnel_Phase_And_Weight(
    const rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const size_t center,
    const size_t offset,
    const double phi_guess,
    double * TMPL_RESTRICT const weight,
    double * TMPL_RESTRICT const psi
)
{
    tmpl_CylFresnelGeometryDouble geo;
    double phi_s;
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
        tau->rho_km_vals[val1], phi_guess
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
    phi_s = tmpl_Double_Stationary_Cyl_Fresnel_Phase_And_Weight(
        tau->k_vals[val0],
        &geo,
        tau->EPS,
        tau->toler,
        weight,
        psi
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
        *psi = -*psi;

    return 57.29577951308232 * phi_s;
}
