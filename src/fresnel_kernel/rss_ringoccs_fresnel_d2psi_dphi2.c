#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double
rssringoccs_Double_Fresnel_d2Psi_dPhi2(double k, double r, double r0,
                                       double phi, double phi0, double B,
                                       double D)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double xi, eta, psi0, dxi, dxi2, deta, deta2, psi_d2, rcpr_D;
    double rcpr_D_squared, cos_B, sin_phi, cos_phi, xi_factor, eta_factor;
    double cos_phi0, sin_phi0, cos_phi_phi0, sin_phi_phi0;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0 / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B   = cos(B);

    cos_phi = cos(phi);
    sin_phi = sin(phi);

    cos_phi0 = cos(phi0);
    sin_phi0 = sin(phi0);

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  This term appears in dxi and dxi2 and xi.                             */
    xi_factor = cos_B * rcpr_D;

    /*  And this term appears in eta, deta, and deta2.                        */
    eta_factor = 2.0 * r * r0 * rcpr_D_squared;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = xi_factor * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - eta_factor*cos_phi_phi0;
    psi0 = sqrt(1.0 + eta - 2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi   = -xi_factor * (r * sin_phi);
    dxi2  = -xi_factor * (r * cos_phi);
    deta  = eta_factor * sin_phi_phi0;
    deta2 = eta_factor * cos_phi_phi0;

    /*  Compute the second partial derivative.                                */
    psi_d2  = (-0.25/(psi0*psi0*psi0)) * (deta - 2.0*dxi) * (deta - 2.0*dxi);
    psi_d2 += (0.5/psi0) * (deta2 - 2.0*dxi2) + dxi2;
    psi_d2 *= k*D;

    return psi_d2;
}

