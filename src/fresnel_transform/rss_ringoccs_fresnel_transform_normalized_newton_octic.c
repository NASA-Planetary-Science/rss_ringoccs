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

#define C00 (+6.720E+02)
#define C01 (-1.680E+02)
#define C02 (+3.200E+01)
#define C03 (-3.000E+00)

#define C10 (+8.064E+03)
#define C11 (-1.008E+03)
#define C12 (+1.280E+02)
#define C13 (-9.000E+00)

#define C20 (-4.880E+02)
#define C21 (+3.380E+02)
#define C22 (-7.200E+01)
#define C23 (+7.000E+00)

#define C30 (-1.952E+03)
#define C31 (+6.760E+02)
#define C32 (-9.600E+01)
#define C33 (+7.000E+00)

#define C40 (+2.900E+01)
#define C41 (-2.600E+01)
#define C42 (+9.000E+00)
#define C43 (-1.000E+00)

#define C50 (+1.160E+02)
#define C51 (-5.200E+01)
#define C52 (+1.200E+01)
#define C53 (-1.000E+00)

#define C60 (-1.400E+01)
#define C61 (+1.400E+01)
#define C62 (-6.000E+00)
#define C63 (+1.000E+00)

#define C70 (-5.600E+01)
#define C71 (+2.800E+01)
#define C72 (-8.000E+00)
#define C73 (+1.000E+00)

#define S0 (1.0 / 105.0)
#define S1 (4.0 / 315.0)
#define S2 (16.0 / 45.0)
#define S3 (64.0 / 45.0)
#define S4 (2048.0 / 45.0)
#define S5 (8192.0 / 45.0)
#define S6 (65536.0 / 315.0)
#define S7 (262144.0 / 315.0)

#define COEFF_POLY_EVAL(z, N) \
(C##N##0*z[0] + C##N##1*z[1] +C##N##2*z[2] +C##N##3*z[3])*S##N

#define ODD_POLY_EVAL(z)    \
z * (coeffs[0] + z##sq*(coeffs[2] + z##sq*(coeffs[4] + z##sq*coeffs[6])))

#define EVEN_POLY_EVAL(z)   \
z##sq * (coeffs[1] + z##sq*(coeffs[3] + z##sq*(coeffs[5] + z##sq*coeffs[7])))

void
rssringoccs_Fresnel_Transform_Normalized_Newton_Octic(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    size_t n, ind[8];

    double scale_factor;
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;
    tmpl_ComplexDouble T_left, T_right, integrand;

    double coeffs[8], diff[4], mean[4], psi[8];

    const size_t shift = nw_pts >> 2;
    const size_t l_ind = center - nw_pts;
    const size_t r_ind = center + nw_pts;

    const double width_actual = 8.0 * tau->dx_km * TMPL_CAST(shift, double);
    const double rcpr_width_actual = 1.0 / width_actual;

    tmpl_ComplexDouble norm = tmpl_CDouble_One;
    tau->T_out[center] = tau->T_in[center];

    ind[0] = center - 4 * shift;
    ind[1] = center - 3 * shift;
    ind[2] = center - 2 * shift;
    ind[3] = center - shift;
    ind[4] = center + shift;
    ind[5] = center + 2 * shift;
    ind[6] = center + 3 * shift;
    ind[7] = center + 4 * shift;


    for (n = 0; n < 8; ++n)
    {
        const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polard(
            tau->rho_km_vals[ind[n]], tau->phi_deg_vals[ind[n]]
        );

        const tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polard(
            tau->rho_km_vals[center], tau->phi_deg_vals[ind[n]]
        );

        const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
            tau->rx_km_vals[ind[n]],
            tau->ry_km_vals[ind[n]],
            tau->rz_km_vals[ind[n]]
        );

        psi[n] = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
        );
    }

    diff[0] = psi[4] - psi[3];
    diff[1] = psi[5] - psi[2];
    diff[2] = psi[6] - psi[1];
    diff[3] = psi[7] - psi[0];

    mean[0] = (psi[4] + psi[3]) * 0.5;
    mean[1] = (psi[5] + psi[2]) * 0.5;
    mean[2] = (psi[6] + psi[1]) * 0.5;
    mean[3] = (psi[7] + psi[0]) * 0.5;

    coeffs[0] = COEFF_POLY_EVAL(diff, 0);
    coeffs[1] = COEFF_POLY_EVAL(mean, 1);
    coeffs[2] = COEFF_POLY_EVAL(diff, 2);
    coeffs[3] = COEFF_POLY_EVAL(mean, 3);
    coeffs[4] = COEFF_POLY_EVAL(diff, 4);
    coeffs[5] = COEFF_POLY_EVAL(mean, 5);
    coeffs[6] = COEFF_POLY_EVAL(diff, 6);
    coeffs[7] = COEFF_POLY_EVAL(mean, 7);

    for (n = 0; n < nw_pts; ++n)
    {
        const double x = x_arr[n] * rcpr_width_actual;
        const double xsq = x * x;

        const double psi_odd = ODD_POLY_EVAL(x);
        const double psi_even = EVEN_POLY_EVAL(x);
        const double psi_left = psi_even + psi_odd;
        const double psi_right = psi_even - psi_odd;

        w_exp_minus_psi_left = tmpl_CDouble_Polar(w_func[n], -psi_left);
        w_exp_minus_psi_right = tmpl_CDouble_Polar(w_func[n], -psi_right);

        T_left = tau->T_in[l_ind + n];
        T_right = tau->T_in[r_ind - n];

        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_left);
        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_right);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_left, T_left);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_right, T_right);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
    }

    scale_factor = tmpl_Double_Rcpr_Sqrt_Two / tmpl_CDouble_Abs(norm);
    integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
