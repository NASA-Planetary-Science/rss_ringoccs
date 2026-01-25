#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>

/*  Numerical integration tools found here.                                   */
#include <libtmpl/include/tmpl_integration.h>

#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

static void rssringoccs_deg4_interp(double * const out, const double * const in)
{
    double mean[2], diff[2];

    mean[0] = (in[1] + in[3]) * 0.5;
    mean[1] = (in[0] + in[4]) * 0.5;
    diff[0] = in[3] - in[1];
    diff[1] = in[4] - in[0];

    out[0] = in[2];
    out[1] = (8.0 * diff[0] - diff[1]) * (1.0 / 3.0);
    out[2] = (16.0 * mean[0] - mean[1]) * (4.0 / 3.0);
    out[3] = (diff[1] - 2.0 * diff[0]) * (16.0 / 3.0);
    out[4] = (mean[1] - 4.0 * mean[0]) * (64.0 / 3.0);
}

void
rssringoccs_Fresnel_Transform_Newton4(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    size_t n;

    double x, xsq;
    tmpl_ComplexDouble l_current, l_next;
    tmpl_ComplexDouble r_current, r_next;
    tmpl_ComplexDouble integrand, in_center;

    double psi_coeffs[5], psi[5], w_coeffs[5], weight[5];
    double psi_odd, psi_even;
    double psi_l_current, psi_l_next;
    double psi_r_current, psi_r_next;
    double w_odd, w_even;
    double w_l_current, w_l_next;
    double w_r_current, w_r_next;

    const size_t shift = nw_pts >> 1;
    const size_t l_ind = center - nw_pts;
    const size_t r_ind = center + nw_pts;
    size_t offset = center - 2 * shift;

    const double width_actual = 4.0 * tau->dx_km * TMPL_CAST(shift, double);
    const double rcpr_width_actual = 1.0 / width_actual;

    tau->T_out[center] = tmpl_CDouble_Zero;

    for (n = 0; n < 5; ++n)
    {
        rssringoccs_Fresnel_Phase_And_Weight(
            tau, center, offset, &weight[n], &psi[n]
        );

        offset += shift;
    }

    rssringoccs_deg4_interp(psi_coeffs, psi);
    rssringoccs_deg4_interp(w_coeffs, weight);

    in_center = tmpl_CDouble_Multiply_Real(weight[2], tau->T_in[center]);

    x = x_arr[0] * rcpr_width_actual;
    xsq = x * x;

    psi_odd = x * (psi_coeffs[1] + xsq * psi_coeffs[3]);
    psi_even = psi_coeffs[0] + xsq * (psi_coeffs[2] + xsq * psi_coeffs[4]);
    psi_l_current = psi_even + psi_odd;
    psi_r_current = psi_even - psi_odd;

    w_odd = x * (w_coeffs[1] + xsq * w_coeffs[3]);
    w_even = w_coeffs[0] + xsq * (w_coeffs[2] + xsq * w_coeffs[4]);
    w_l_current = (w_even + w_odd) * w_func[0];
    w_r_current = (w_even - w_odd) * w_func[0];

    l_current = tmpl_CDouble_Multiply_Real(w_l_current, tau->T_in[l_ind]);
    r_current = tmpl_CDouble_Multiply_Real(w_r_current, tau->T_in[r_ind]);

    for (n = 1; n < nw_pts; ++n)
    {
        x = x_arr[n] * rcpr_width_actual;
        xsq = x * x;

        psi_odd = x * (psi_coeffs[1] + xsq * psi_coeffs[3]);
        psi_even = psi_coeffs[0] + xsq * (psi_coeffs[2] + xsq * psi_coeffs[4]);
        psi_l_next = psi_even + psi_odd;
        psi_r_next = psi_even - psi_odd;

        w_odd = x * (w_coeffs[1] + xsq * w_coeffs[3]);
        w_even = w_coeffs[0] + xsq * (w_coeffs[2] + xsq * w_coeffs[4]);
        w_l_next = (w_even + w_odd) * w_func[n];
        w_r_next = (w_even - w_odd) * w_func[n];

        l_next = tmpl_CDouble_Multiply_Real(w_l_next, tau->T_in[l_ind + n]);
        r_next = tmpl_CDouble_Multiply_Real(w_r_next, tau->T_in[r_ind - n]);

        integrand = tmpl_CDouble_Filon11_Integrand(
            l_current, l_next, psi_l_current, psi_l_next, tau->dx_km
        );

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Filon11_Integrand(
            r_next, r_current, psi_r_next, psi_r_current, tau->dx_km
        );

        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        l_current = l_next;
        r_current = r_next;
        psi_l_current = psi_l_next;
        psi_r_current = psi_r_next;
        w_l_current = w_l_next;
        w_r_current = w_r_next;
    }

    integrand = tmpl_CDouble_Filon12_Integrand(
        l_current,
        in_center,
        r_current,
        psi_l_current,
        psi[2],
        psi_r_current,
        tau->dx_km
    );

    tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
}
