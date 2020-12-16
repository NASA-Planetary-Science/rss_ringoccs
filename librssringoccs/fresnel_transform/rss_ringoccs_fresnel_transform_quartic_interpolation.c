#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stdlib.h>

rssringoccs_ComplexDouble
Fresnel_Transform_Quartic_Double(double *x_arr, double *phi_arr,
                                 rssringoccs_ComplexDouble *T_in,
                                 double *w_func, double kD, double r,
                                 double B, double D,
                                 double EPS, unsigned long toler, double dx,
                                 double F, unsigned long n_pts,
                                 unsigned long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3], rcpr_F, rcpr_w, rcpr_w_sq;
    double* psi_n;
    double psi, phi, width, sin_psi, cos_psi;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;

    psi_n = (double *)malloc(sizeof(double)*n_pts);
    width = x_arr[n_pts-1] - x_arr[0];
    rcpr_w = 1.0 / width;
    rcpr_w_sq = rcpr_w * rcpr_w;
    rcpr_F = 1.0/F;


    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;

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
    psi_half_diff = psi_n[(n_pts-1)/4] - psi_n[3*(n_pts-1)/4];
    psi_full_diff = psi_n[0] - psi_n[n_pts-1];

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[1] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;
    C[2] = (psi_full_mean-4.0*psi_half_mean)*rcpr_w_sq*rcpr_w_sq*21.33333333333;

    for (i = 0; i<n_pts; ++i){

        psi = C[2];
        psi = psi*(r-x_arr[i])+C[1];
        psi = psi*(r-x_arr[i])+C[0];
        psi = psi*(r-x_arr[i])*(r-x_arr[i]);

        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[i], exp_psi);
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[j]);
        T_out = rssringoccs_CDouble_Add(T_out, integrand);
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(0.5*dx*rcpr_F, 0.5*dx*rcpr_F);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
