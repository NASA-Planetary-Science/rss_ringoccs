#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>

#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

#define ONE_THIRD (3.333333333333333333333333333333E-01)
#define FOUR_THIRDS (1.333333333333333333333333333333E+00)
#define SIXTEEN_THIRDS (5.333333333333333333333333333333E+00)
#define SIXTY_FOUR_THIRDS (2.133333333333333333333333333333E+01)

void
rssringoccs_Fresnel_Transform_Newton_Quartic_Norm(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const w_func,
    size_t n_pts,
    size_t center
)
{
    /*  Variable for indexing.                                                */
    size_t n, ind[4];

    double x, xsq, scale_factor, width, rcpr_width;
    double psi_odd, psi_even;
    double psi_left, psi_right;
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;
    tmpl_ComplexDouble T_left, T_right, integrand;
    tmpl_ComplexDouble norm = tmpl_CDouble_Zero;

    /*  Cast constants to size_t to avoid implicit conversion.                */
    const size_t zero = TMPL_CAST(0, size_t);
    const size_t one = TMPL_CAST(1, size_t);
    const size_t four = TMPL_CAST(4, size_t);

    double coeffs[4];
    double psi[4], half_diff, full_diff, half_mean, full_mean;

    const size_t half_n_pts = (n_pts - one) >> 1;
    const size_t shift = half_n_pts >> 1;
    const size_t l_ind = center - half_n_pts;
    const size_t r_ind = center + half_n_pts;

    ind[0] = center - 2 * shift;
    ind[1] = center - shift;
    ind[2] = center + shift;
    ind[3] = center + 2 * shift;

    width = tau->rho_km_vals[ind[3]] - tau->rho_km_vals[ind[0]];
    rcpr_width = 1.0 / width;

    tau->T_out[center] = tmpl_CDouble_Zero;

    for (n = zero; n < four; ++n)
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

    half_mean = (psi[1] + psi[2]) * 0.5;
    full_mean = (psi[0] + psi[3]) * 0.5;
    half_diff = psi[1] - psi[2];
    full_diff = psi[0] - psi[3];

    coeffs[0] = (8.0 * half_diff - full_diff) * ONE_THIRD;
    coeffs[1] = (16.0 * half_mean - full_mean) * FOUR_THIRDS;
    coeffs[2] = (full_diff - 2.0 * half_diff) * SIXTEEN_THIRDS;
    coeffs[3] = (full_mean - 4.0 * half_mean) * SIXTY_FOUR_THIRDS;

    for (n = zero; n < half_n_pts; ++n)
    {
        x = (tau->rho_km_vals[center] - tau->rho_km_vals[l_ind + n])*rcpr_width;
        xsq = x * x;

        psi_odd = x * (coeffs[0] + xsq * coeffs[2]);
        psi_even = xsq * (coeffs[1] + xsq * coeffs[3]);
        psi_left = psi_even + psi_odd;
        psi_right = psi_even - psi_odd;

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

    /*  Add the central point in the Riemann sum. This is the center of the   *
     *  window function. That is, where w_func = 1.                           */
    tmpl_CDouble_AddTo(&tau->T_out[center], &tau->T_in[center]);
    tmpl_CDouble_AddTo_Real(&norm, 1.0);

    /*  The integral in the numerator of the scale factor is F sqrt(2),       *
     *  the denominator is |norm| dx, so the ratio is sqrt(2) F / |norm| dx.  *
     *  Since the dx is on the bottom, this will cancel the dx that occurs in *
     *  the Riemann sum for T_out. The complex scale factor for the Fresnel   *
     *  transform (the constant term outside the integral) is (1 + i) / 2 F.  *
     *  We therefore have:                                                    *
     *                                                                        *
     *      1 + i numer      1 + i sqrt(2) F                                  *
     *      ----- ----- dx = ----- --------- dx                               *
     *       2 F  demom       2 F  |norm| dx                                  *
     *                                                                        *
     *                            1 + i                                       *
     *                     = --------------                                   *
     *                       sqrt(2) |norm|                                   *
     *                                                                        *
     *  Compute this and scale the result to finish the calculation.          */
    scale_factor = tmpl_Double_Rcpr_Sqrt_Two / tmpl_CDouble_Abs(norm);
    integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}

#undef ONE_THIRD
#undef FOUR_THIRDS
#undef SIXTEEN_THIRDS
#undef SIXTY_FOUR_THIRDS
