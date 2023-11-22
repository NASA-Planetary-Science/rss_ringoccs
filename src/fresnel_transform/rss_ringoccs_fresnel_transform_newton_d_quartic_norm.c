#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

#define ONE_THIRD (0.3333333333333333333333333333)
#define FOUR_THIRDS (1.33333333333333333333333333)
#define SIXTEEN_THIRDS (5.3333333333333333333333333)
#define SIXTY_FOUR_THIRDS (21.3333333333333333333333333)

void
rssringoccs_Fresnel_Transform_Newton_D_Quartic_Norm(rssringoccs_TAUObj *tau,
                                                    const double *w_func,
                                                    size_t n_pts,
                                                    size_t center)
{
    /*  Variable for indexing.                                                */
    size_t n;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[4], abs_norm, real_norm;
    double psi_n[4], psi_half_diff, psi_full_diff, psi, phi;
    double psi_half_mean, psi_full_mean, D, x;
    tmpl_ComplexDouble exp_psi, norm, integrand;

    const double rcpr_w = 1.0 / tau->w_km_vals[center];
    const double rcpr_w_sq = rcpr_w * rcpr_w;
    const double rcpr_w_cb = rcpr_w_sq * rcpr_w;
    const double rcpr_w_qr = rcpr_w_sq * rcpr_w_sq;
    size_t offset = center - ((n_pts-1) >> 1U);

    const double rho[4] = {
        tau->rho_km_vals[center] - 0.5*tau->w_km_vals[center],
        tau->rho_km_vals[center] - 0.25*tau->w_km_vals[center],
        tau->rho_km_vals[center] + 0.25*tau->w_km_vals[center],
        tau->rho_km_vals[center] + 0.5*tau->w_km_vals[center]
    };

    const double num = tau->phi_deg_vals[offset + n_pts - 1] - tau->phi_deg_vals[offset];
    const double den = tau->rho_km_vals[offset + n_pts - 1] - tau->rho_km_vals[offset];
    const double slope = num / den;

    const double phi0[4] = {
        (rho[0] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset],
        (rho[1] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset],
        (rho[2] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset],
        (rho[3] - tau->rho_km_vals[offset])*slope + tau->phi_deg_vals[offset]
    };

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    tau->T_out[center] = tmpl_CDouble_Zero;
    norm = tmpl_CDouble_Zero;

    for (n = (size_t)0; n < (size_t)4; ++n)
    {
        phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_D_Newton(
            tau->k_vals[center],        /*  Wavenumber.                       */
            tau->rho_km_vals[center],   /*  Ring radius.                      */
            rho[n],                     /*  Dummy radius.                     */
            phi0[n],                    /*  Initial stationary azimuth guess. */
            phi0[n],                    /*  Dummy azimuthal angle.            */
            tau->B_deg_vals[center],    /*  Ring opening angle.               */
            tau->rx_km_vals[center],    /*  x-coordinate of spacecraft.       */
            tau->ry_km_vals[center],    /*  y-coordinate of spacecraft.       */
            tau->rz_km_vals[center],    /*  z-coordinate of spacecraft.       */
            tau->EPS,                   /*  Epsilon error for Newton's method.*/
            tau->toler                  /*  Maximum number of iterations.     */
        );

        D = tmpl_Double_Cyl_Fresnel_Observer_Distance(
            rho[n],                     /* Ring radius.                       */
            phi,                        /* Stationary azimuth angle.          */
            tau->rx_km_vals[center],    /* Cassini x coordinate.              */
            tau->ry_km_vals[center],    /* Cassini y coordinate.              */
            tau->rz_km_vals[center]     /* Cassini z coordinate.              */
        );

        psi_n[n] = tmpl_Double_Cyl_Fresnel_Psi(
            tau->k_vals[center],        /*  Wavenumber.                       */
            tau->rho_km_vals[center],   /*  Ring radius.                      */
            rho[n],                     /*  Dummy radius.                     */
            phi,                        /*  Stationary azimuthal angle.       */
            phi0[n],                    /*  Dummy azimuthal angle.            */
            tau->B_deg_vals[center],    /*  Ring opening angle.               */
            D                           /*  Ring-Spacecraft distance.         */
        );
    }

    psi_half_mean = (psi_n[1] + psi_n[2]) * 0.5;
    psi_full_mean = (psi_n[0] + psi_n[3]) * 0.5;
    psi_half_diff = psi_n[1] - psi_n[2];
    psi_full_diff = psi_n[0] - psi_n[3];

    C[0] = (8.0*psi_half_diff - psi_full_diff) * rcpr_w * ONE_THIRD;
    C[1] = (16.0*psi_half_mean - psi_full_mean) * rcpr_w_sq * FOUR_THIRDS;
    C[2] = (psi_full_diff - 2.0*psi_half_diff) * rcpr_w_cb * SIXTEEN_THIRDS;
    C[3] = (psi_full_mean - 4.0*psi_half_mean) * rcpr_w_qr * SIXTY_FOUR_THIRDS;

    for (n = (size_t)0; n < n_pts; ++n)
    {
        x = tau->rho_km_vals[center] - tau->rho_km_vals[offset];
        psi = x*(C[0] + x*(C[1] + x*(C[2] + x*C[3])));

        exp_psi = tmpl_CDouble_Polar(w_func[n], -psi);
        integrand = tmpl_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        tmpl_CDouble_AddTo(&norm, &exp_psi);
        offset++;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    abs_norm = tmpl_CDouble_Abs(norm);
    real_norm = tmpl_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = tmpl_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
