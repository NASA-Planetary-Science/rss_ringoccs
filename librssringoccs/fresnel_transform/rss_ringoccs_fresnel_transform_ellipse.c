#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

rssringoccs_ComplexDouble
Fresnel_Transform_Ellipse_Double(rssringoccs_TAUObj *tau,
                                 double *w_func,
                                 unsigned long n_pts,
                                 unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long m, offset;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi, cos_psi, sin_psi, factor, x, y, z, dx, dy, D;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    factor = 0.5 * tau->dx_km / tau->F_km_vals[center];

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    offset = center-(unsigned long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m < n_pts; ++m)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = rssringoccs_Double_Newton_Raphson_Fresnel_Ellipse(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset],
            tau->phi_rad_vals[offset],
            tau->phi_rad_vals[offset],
            tau->B_rad_vals[center],
            tau->ecc,
            tau->peri,
            tau->EPS,
            tau->toler,
            tau->rx_km_vals[center],
            tau->ry_km_vals[center],
            tau->rz_km_vals[center]
        );

        x = tau->rho_km_vals[offset] * cos(phi);
        y = tau->rho_km_vals[offset] * sin(phi);
        z = tau->rz_km_vals[center];
        dx = x - tau->rx_km_vals[center];
        dy = y - tau->ry_km_vals[center];
        D = sqrt(dx*dx + dy*dy + z*z);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = rssringoccs_Double_Fresnel_Psi(
            tau->k_vals[center],
            tau->rho_km_vals[center],
            tau->rho_km_vals[offset],
            phi,
            tau->phi_rad_vals[offset],
            tau->B_rad_vals[center],
            D
        );

        cos_psi = w_func[m]*rssringoccs_Double_Cos(psi);
        sin_psi = w_func[m]*rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, tau->T_in[offset]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        offset += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(factor, factor);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
