#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

rssringoccs_ComplexDouble
Fresnel_Transform_Ellipse_Double(double *x_arr, double *phi_arr,
                                 rssringoccs_ComplexDouble *T_in,
                                 double *w_func, double kD, double r, double B,
                                 double D, double EPS, unsigned long toler,
                                 double dx, double F, unsigned long n_pts,
                                 unsigned long center, double ecc, double peri)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi, cos_psi, sin_psi, rcpr_F;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    rcpr_F = 1.0/F;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = rssringoccs_Double_Newton_Raphson_Fresnel_Ellipse(
            kD, r, x_arr[i], phi_arr[i], phi_arr[i], B, D, ecc, peri, EPS, toler
        );

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi = rssringoccs_Double_Fresnel_Psi(kD, r, x_arr[i], phi, phi_arr[i], B, D);

        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[i], exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[i]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(0.5*dx*rcpr_F, 0.5*dx*rcpr_F);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
