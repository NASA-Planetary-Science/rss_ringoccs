/*  TMPL_RESTRICT macro provided here.                                        */
#include <libtmpl/include/tmpl_config.h>

/*  Double precision complex numbers and routines given here.                 */
#include <libtmpl/include/tmpl_complex.h>

/*  Numerical integration tools found here.                                   */
#include <libtmpl/include/tmpl_integration.h>

/*  rssringoccs_TAUObj typedef provided here.                                 */
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>

/*  Function prototype / forward declaration found here.                      */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Inverse Fresnel transform via Newton-Raphson with window normalization.   */
void
rssringoccs_Fresnel_Transform_Newton_Filon12(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    const size_t nw_pts,
    const size_t center
)
{
    double weight, left_scale, mid_scale, right_scale;
    double left_psi, mid_psi, right_psi;

    tmpl_ComplexDouble left, mid, right, integrand;

    size_t offset = center - (nw_pts >> 1);
    size_t n;

    tau->T_out[center] = tmpl_CDouble_Zero;

    rssringoccs_Fresnel_Phase_And_Weight(
        tau, center, offset, &weight, &left_psi
    );

    left_scale = weight * w_func[0];
    left = tmpl_CDouble_Multiply_Real(left_scale, tau->T_in[offset]);

    for (n = 0; n < nw_pts - 2; n += 2)
    {
        rssringoccs_Fresnel_Phase_And_Weight(
            tau, center, offset + 1, &weight, &mid_psi
        );

        mid_scale = weight * w_func[n + 1];
        mid = tmpl_CDouble_Multiply_Real(mid_scale, tau->T_in[offset + 1]);

        rssringoccs_Fresnel_Phase_And_Weight(
            tau, center, offset + 2, &weight, &right_psi
        );

        right_scale = weight * w_func[n + 2];
        right = tmpl_CDouble_Multiply_Real(right_scale, tau->T_in[offset + 2]);

        integrand = tmpl_CDouble_Filon12_Integrand(
            left,
            mid,
            right,
            left_psi,
            mid_psi,
            right_psi,
            tau->dx_km
        );

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        left_scale = right_scale;
        left_psi = right_psi;
        left = right;
        offset += 2;
    }
}
/*  End of rssringoccs_Fresnel_Transform_Newton_Filon12.                      */
