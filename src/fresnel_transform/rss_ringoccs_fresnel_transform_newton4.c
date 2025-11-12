#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

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
    size_t ind[4];

    double scale_factor;
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;
    tmpl_ComplexDouble T_left, T_right, integrand;
    tmpl_CylFresnelGeometryDouble geo;

    double coeffs[4];
    double psi[4], diff[2], mean[2];

    const size_t shift = nw_pts >> 1;
    const size_t l_ind = center - nw_pts;
    const size_t r_ind = center + nw_pts;

    const double width_actual = 4.0 * tau->dx_km * TMPL_CAST(shift, double);
    const double rcpr_width_actual = 1.0 / width_actual;

    tmpl_ComplexDouble norm = tmpl_CDouble_One;
    tau->T_out[center] = tau->T_in[center];

    ind[0] = center - 2 * shift;
    ind[1] = center - shift;
    ind[2] = center + shift;
    ind[3] = center + 2 * shift;

    for (n = 0; n < 4; ++n)
    {
        geo.position = tmpl_3DDouble_Rect(
            tau->rx_km_vals[ind[n]],
            tau->ry_km_vals[ind[n]],
            tau->rz_km_vals[ind[n]]
        );

        geo.intercept = tmpl_2DDouble_Polard(
            tau->rho_km_vals[ind[n]], tau->phi_deg_vals[ind[n]]
        );

        geo.dummy = tmpl_2DDouble_Polard(
            tau->rho_km_vals[center], tau->phi_deg_vals[ind[n]]
        );

        psi[n] = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &geo, tau->EPS, tau->toler
        );
    }

    mean[0] = (psi[1] + psi[2]) * 0.5;
    mean[1] = (psi[0] + psi[3]) * 0.5;
    diff[0] = psi[2] - psi[1];
    diff[1] = psi[3] - psi[0];

    coeffs[0] = (8.0 * diff[0] - diff[1]) * (1.0 / 3.0);
    coeffs[1] = (16.0 * mean[0] - mean[1]) * (4.0 / 3.0);
    coeffs[2] = (diff[1] - 2.0 * diff[0]) * (16.0 / 3.0);
    coeffs[3] = (mean[1] - 4.0 * mean[0]) * (64.0 / 3.0);

    for (n = 0; n < nw_pts; ++n)
    {
        const double x = x_arr[n] * rcpr_width_actual;
        const double xsq = x * x;

        const double psi_odd = x * (coeffs[0] + xsq * coeffs[2]);
        const double psi_even = xsq * (coeffs[1] + xsq * coeffs[3]);
        const double psi_left = psi_even + psi_odd;
        const double psi_right = psi_even - psi_odd;

        w_exp_minus_psi_left = tmpl_CDouble_Polar(w_func[n], -psi_left);
        w_exp_minus_psi_right = tmpl_CDouble_Polar(w_func[n], -psi_right);

        T_left = tau->T_in[l_ind + n];
        T_right = tau->T_in[r_ind - n];

        if (tau->use_norm)
        {
            tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_left);
            tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_right);
        }

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_left, T_left);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_right, T_right);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
    }

    if (tau->use_norm)
    {
        scale_factor = tmpl_Double_Rcpr_Sqrt_Two / tmpl_CDouble_Abs(norm);
        integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
        tmpl_CDouble_MultiplyBy(&tau->T_out[center], &integrand);
    }
}
