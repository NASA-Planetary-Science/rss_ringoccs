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

#define C00 (+6.400000000000000000000000E+0)
#define C01 (-1.600000000000000000000000E+0)
#define C02 (+3.047619047619047619047620E-1)
#define C03 (-2.857142857142857142857140E-2)

#define C10 (+1.024000000000000000000000E+2)
#define C11 (-1.280000000000000000000000E+1)
#define C12 (+1.625396825396825396825400E+0)
#define C13 (-1.142857142857142857142860E-1)

#define C20 (-1.735111111111111111111110E+2)
#define C21 (+1.201777777777777777777780E+2)
#define C22 (-2.560000000000000000000000E+1)
#define C23 (+2.488888888888888888888890E+0)

#define C30 (-2.776177777777777777777780E+3)
#define C31 (+9.614222222222222222222220E+2)
#define C32 (-1.365333333333333333333330E+2)
#define C33 (+9.955555555555555555555560E+0)

#define C40 (+1.319822222222222222222220E+3)
#define C41 (-1.183288888888888888888890E+3)
#define C42 (+4.096000000000000000000000E+2)
#define C43 (-4.551111111111111111111110E+1)

#define C50 (+2.111715555555555555555560E+4)
#define C51 (-9.466311111111111111111110E+3)
#define C52 (+2.184533333333333333333330E+3)
#define C53 (-1.820444444444444444444440E+2)

#define C60 (-2.912711111111111111111110E+3)
#define C61 (+2.912711111111111111111110E+3)
#define C62 (-1.248304761904761904761900E+3)
#define C63 (+2.080507936507936507936510E+2)

#define C70 (-4.660337777777777777777780E+4)
#define C71 (+2.330168888888888888888890E+4)
#define C72 (-6.657625396825396825396830E+3)
#define C73 (+8.322031746031746031746030E+2)

#define COEFF_POLY_EVAL(z, N) \
(C##N##0*z[0] + C##N##1*z[1] +C##N##2*z[2] +C##N##3*z[3])

#define ODD_POLY_EVAL(z)    \
z * (coeffs[0] + z##sq*(coeffs[2] + z##sq*(coeffs[4] + z##sq*coeffs[6])))

#define EVEN_POLY_EVAL(z)   \
z##sq * (coeffs[1] + z##sq*(coeffs[3] + z##sq*(coeffs[5] + z##sq*coeffs[7])))

void
rssringoccs_Fresnel_Transform_Newton8(
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
    tmpl_CylFresnelGeometryDouble geo;

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
