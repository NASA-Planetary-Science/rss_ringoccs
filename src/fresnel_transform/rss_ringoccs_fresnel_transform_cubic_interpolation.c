#include <math.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

void
rssringoccs_Fresnel_Transform_Cubic(rssringoccs_TAUObj *tau,
                                    double *w_func,
                                    unsigned long n_pts,
                                    unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, ind[4], offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3], factor, rcpr_w, rcpr_w_sq, psi_n[4], x;
    double psi, phi, sin_psi, cos_psi;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    tmpl_ComplexDouble exp_psi, integrand;

    rcpr_w = 1.0 / tau->w_km_vals[center];
    rcpr_w_sq = rcpr_w * rcpr_w;
    factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    ind[0] = 0;
    ind[1] = (n_pts-1)/4;
    ind[2] = 3*(n_pts-1)/4;
    ind[3] = n_pts - 1;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center-(unsigned long)((n_pts-1)/2);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.       */
    for (i = 0; i < 4; ++i)
    {

        phi = rssringoccs_Newton_Raphson_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            tau->D_km_vals[center],
            tau->EPS,
            tau->toler
        );

        psi_n[i] = rssringoccs_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset + ind[i]],
            phi,
            tau->phi_rad_vals[offset + ind[i]],
            tau->B_rad_vals[center],
            tau->D_km_vals[center]
        );
    }

    psi_half_mean = (psi_n[1] + psi_n[2]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[3]) / 2;
    psi_half_diff = psi_n[1] - psi_n[2];
    psi_full_diff = psi_n[0] - psi_n[3];

    C[0] = (8*psi_half_diff - psi_full_diff) * rcpr_w * 0.333333333333333333333;
    C[1] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[2] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;

    for (i = 0; i<n_pts; ++i)
    {
        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];
        psi = C[2];
        psi = psi*x + C[1];
        psi = psi*x + C[0];
        psi = psi*x;

        cos_psi = w_func[i]*cos(psi);
        sin_psi = w_func[i]*sin(psi);
        exp_psi = tmpl_CDouble_Rect(cos_psi, -sin_psi);
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tau->T_out[center] = tmpl_CDouble_Add(tau->T_out[center], integrand);
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand,
                                                      tau->T_out[center]);
}
