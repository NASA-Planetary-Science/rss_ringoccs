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
rssringoccs_Fresnel_Transform_Newton_Filon01(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    const size_t nw_pts,
    const size_t center
)
{
    double weight, left_scale, right_scale;
    double left_psi, right_psi;

    tmpl_ComplexDouble left, right, integrand, midpoint;
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;
    size_t offset = center - (nw_pts >> 1);
    size_t n;

    tau->T_out[center] = tmpl_CDouble_Zero;

    rssringoccs_Fresnel_Phase_And_Weight(
        tau, center, offset, &weight, &left_psi
    );

    left_scale = weight * w_func[0];
    left = tmpl_CDouble_Multiply_Real(left_scale, tau->T_in[offset]);

    for (n = 0; n < nw_pts - 1; ++n)
    {
        rssringoccs_Fresnel_Phase_And_Weight(
            tau, center, offset + 1, &weight, &right_psi
        );

        right_scale = weight * w_func[n + 1];
        right = tmpl_CDouble_Multiply_Real(right_scale, tau->T_in[offset + 1]);

        midpoint = tmpl_CDouble_Midpoint(left, right);

        integrand = tmpl_CDouble_Filon01_Integrand(
            midpoint, left_psi, right_psi, tau->dx_km
        );

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        if (tau->use_norm)
        {
            const double average = 0.5 * (left_scale + right_scale);
            integrand = tmpl_Double_Filon01_Integrand(
                average, left_psi, right_psi, tau->dx_km
            );

            tmpl_CDouble_AddTo(&norm, &integrand);
        }

        left_scale = right_scale;
        left_psi = right_psi;
        left = right;
        ++offset;
    }

    /*  By Babinet's principle, the free space integral is 1. We are          *
     *  integrating over a finite data set, and have introduced a tapering    *
     *  function. The normalization factor is the magnitude of the free space *
     *  integral across the entire real line divided by the tapered integral  *
     *  across the window, which is hence 1 / | norm |. Compute this.         */
    if (tau->use_norm)
    {
        /*  Scale factor computed using Babinet's principle.                  */
        const double scale = 1.0 / tmpl_CDouble_Abs(norm);

        /*  Scale the Riemann sum by the normalization factor to conclude.    */
        tmpl_CDouble_MultiplyBy_Real(&tau->T_out[center], scale);
    }
}
/*  End of rssringoccs_Fresnel_Transform_Newton_Filon01.                      */
