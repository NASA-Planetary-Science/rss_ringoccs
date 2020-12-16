#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stdlib.h>

rssringoccs_ComplexDouble
Fresnel_Transform_Quadratic_Norm_Double(double *x_arr, double *phi_arr,
                                        rssringoccs_ComplexDouble *T_in,
                                        double *w_func, double kD, double r,
                                        double B, double D, double EPS,
                                        unsigned long toler,
                                        unsigned long n_pts,
                                        unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C_2, abs_norm, real_norm;
    double *psi_n;
    double psi, phi, width, rcpr_w, rcpr_w_sq;
    double psi_half_mean, psi_full_mean, cos_psi, sin_psi;
    rssringoccs_ComplexDouble T_out, exp_psi, norm, integrand;

    /*  Allocate memory for psi_n. Per C99 recommendations, we do not cast    *
     *  malloc since void is safely promoted to whatever type psi_n is.       */
    psi_n = malloc(sizeof(*psi_n) * n_pts);

    width = x_arr[n_pts-1] - x_arr[0];
    rcpr_w = 1.0 / width;
    rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    norm  = rssringoccs_CDouble_Zero;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i)
    {
        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = rssringoccs_Double_Fresnel_Psi(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;
    free(psi_n);

    C_2 = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;

    for (i = 0; i<n_pts; ++i)
    {
        psi = (r-x_arr[i])*(r-x_arr[i])*C_2;

        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[i], exp_psi);
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[j]);
        T_out = rssringoccs_CDouble_Add(T_out, integrand);
        norm = rssringoccs_CDouble_Add(norm, exp_psi);
        j     += 1;
    }

    /*  Add the central point in the Riemann sum. This is center of the       *
     *  window function. That is, where w_func = 1.                           */
    T_out = rssringoccs_CDouble_Add(T_out, T_in[center]);
    norm  = rssringoccs_CDouble_Add_Real(1.0, norm);
    abs_norm = rssringoccs_CDouble_Abs(norm);
    real_norm = rssringoccs_Sqrt_Two / abs_norm;

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    integrand = rssringoccs_CDouble_Rect(0.5*real_norm, 0.5*real_norm);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
