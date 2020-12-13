#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float
rssringoccs_Float_Fresnel_dPsi_dPhi(float kD, float r, float r0, float phi,
                                    float phi0, float B, float D)
{
    float xi, eta, psi0, dxi, deta, cos_B, rcpr_D, rcpr_D_squared, cos_B_by_D;
    float cos_phi, sin_phi, cos_phi0, sin_phi0, cos_phi_phi0, sin_phi_phi0;
    float factor, dpsi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0F / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B   = rssringoccs_Float_Cos(B);

    cos_phi = rssringoccs_Float_Cos(phi);
    sin_phi = rssringoccs_Float_Sin(phi);

    cos_phi0 = rssringoccs_Float_Cos(phi0);
    sin_phi0 = rssringoccs_Float_Sin(phi0);

    /*  These terms occur frequently.                                         */
    cos_B_by_D = cos_B*rcpr_D;
    factor = 2.0F * r * r0 * rcpr_D_squared;

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cos_B_by_D * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - factor*cos_phi_phi0;
    psi0 = rssringoccs_Float_Sqrt(1.0F + eta - 2.0F*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cos_B_by_D * (r * sin_phi);
    deta = factor * sin_phi_phi0;
    dpsi = kD * ((0.5F/psi0)*(deta - 2.0F*dxi) + dxi);

    return dpsi;
}

double
rssringoccs_Double_Fresnel_dPsi_dPhi(double kD, double r, double r0, double phi,
                                     double phi0, double B, double D)
{
    double xi, eta, psi0, dxi, deta, cos_B, rcpr_D, rcpr_D_squared, cos_B_by_D;
    double cos_phi, sin_phi, cos_phi0, sin_phi0, cos_phi_phi0, sin_phi_phi0;
    double factor, dpsi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0 / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B   = rssringoccs_Double_Cos(B);

    cos_phi = rssringoccs_Double_Cos(phi);
    sin_phi = rssringoccs_Double_Sin(phi);

    cos_phi0 = rssringoccs_Double_Cos(phi0);
    sin_phi0 = rssringoccs_Double_Sin(phi0);

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
    psi0 = rssringoccs_Double_Sqrt(1.0 + eta - 2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cos_B_by_D * (r * sin_phi);
    deta = factor * sin_phi_phi0;
    dpsi = kD * ((0.5/psi0)*(deta - 2.0*dxi) + dxi);

    return dpsi;
}

long double
rssringoccs_LDouble_Fresnel_dPsi_dPhi(long double kD, long double r,
                                     long double r0, long double phi,
                                     long double phi0, long double B,
                                     long double D)
{
    long double xi, eta, psi0, dxi, deta, cos_B, rcpr_D, rcpr_D_squared;
    long double cos_B_by_D, cos_phi, sin_phi, cos_phi0, sin_phi0, cos_phi_phi0;
    long double sin_phi_phi0, factor, dpsi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0L / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B   = rssringoccs_LDouble_Cos(B);

    cos_phi = rssringoccs_LDouble_Cos(phi);
    sin_phi = rssringoccs_LDouble_Sin(phi);

    cos_phi0 = rssringoccs_LDouble_Cos(phi0);
    sin_phi0 = rssringoccs_LDouble_Sin(phi0);

    /*  These terms occur frequently.                                         */
    cos_B_by_D = cos_B*rcpr_D;
    factor = 2.0L * r * r0 * rcpr_D_squared;

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cos_B_by_D * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - factor*cos_phi_phi0;
    psi0 = rssringoccs_LDouble_Sqrt(1.0L + eta - 2.0L*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cos_B_by_D * (r * sin_phi);
    deta = factor * sin_phi_phi0;
    dpsi = kD * ((0.5L/psi0)*(deta - 2.0L*dxi) + dxi);

    return dpsi;
}
