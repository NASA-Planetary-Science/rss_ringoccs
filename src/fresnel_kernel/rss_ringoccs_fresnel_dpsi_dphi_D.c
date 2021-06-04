#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double
rssringoccs_Double_Fresnel_dPsi_dPhi_D(double k, double r, double r0, double phi,
                                       double phi0, double B, double rx,
                                       double ry, double rz)
{
    double xi, eta, psi0, dxi, deta, cos_B, rcpr_D, rcpr_D_squared, cos_B_by_D;
    double cos_phi, sin_phi, cos_phi0, sin_phi0, cos_phi_phi0, sin_phi_phi0;
    double factor, dpsi, dD,D, dx, dy;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B   = cos(B);

    cos_phi = cos(phi);
    sin_phi = sin(phi);

    cos_phi0 = cos(phi0);
    sin_phi0 = sin(phi0);

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    dx = rx - r*cos_phi;
    dy = ry - r*sin_phi;
    D = sqrt(dx*dx + dy*dy + rz*rz);
    rcpr_D = 1.0 / D;
    rcpr_D_squared = rcpr_D * rcpr_D;
    dD = r*(dx*sin_phi-dy*cos_phi)*rcpr_D;

    /*  These terms occur frequently.                                         */
    cos_B_by_D = cos_B*rcpr_D;
    factor = 2.0 * r * r0 * rcpr_D_squared;

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cos_B_by_D * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - factor*cos_phi_phi0;
    psi0 = sqrt(1.0 + eta - 2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi   = -cos_B_by_D * (r * sin_phi) - xi*rcpr_D*dD;
    deta  = factor * sin_phi_phi0 - 2.0*eta*rcpr_D*dD;
    dpsi  = k*D * ((0.5/psi0)*(deta - 2.0*dxi) + dxi);
    dpsi += k*dD*(psi0 + xi - 1.0);

    return dpsi;
}
