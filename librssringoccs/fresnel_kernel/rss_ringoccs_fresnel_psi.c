#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float
rssringoccs_Float_Fresnel_Psi(float kD, float r, float r0, float phi,
                              float phi0, float B, float D)
{
    float xi, eta, rcpr_D, rcpr_D_squared, cos_B, cos_phi, cos_phi0;
    float cos_phi_phi0, psi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0F / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B        = rssringoccs_Float_Cos(B);
    cos_phi      = rssringoccs_Float_Cos(phi);
    cos_phi0     = rssringoccs_Float_Cos(phi0);
    cos_phi_phi0 = rssringoccs_Float_Cos(phi-phi0);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos_B * rcpr_D) * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r - 2.0F*r*r0*cos_phi_phi0) * rcpr_D_squared;

    /* Sign of xi swapped from MTR86.                                         */
    psi = kD * (rssringoccs_Float_Sqrt(1.0F + eta - 2.0F*xi) + xi - 1.0F);

    return psi;
}

double
rssringoccs_Double_Fresnel_Psi(double kD, double r, double r0, double phi,
                               double phi0, double B, double D)
{
    double xi, eta, rcpr_D, rcpr_D_squared, cos_B, cos_phi, cos_phi0;
    double cos_phi_phi0, psi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0 / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B        = rssringoccs_Double_Cos(B);
    cos_phi      = rssringoccs_Double_Cos(phi);
    cos_phi0     = rssringoccs_Double_Cos(phi0);
    cos_phi_phi0 = rssringoccs_Double_Cos(phi-phi0);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos_B * rcpr_D) * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos_phi_phi0) * rcpr_D_squared;

    /* Sign of xi swapped from MTR86.                                         */
    psi = kD * (rssringoccs_Double_Sqrt(1.0 + eta - 2.0*xi) + xi - 1.0);

    return psi;
}

long double
rssringoccs_LDouble_Fresnel_Psi(long double kD, long double r, long double r0,
                                long double phi, long double phi0, long double B,
                                long double D)
{
    long double xi, eta, rcpr_D, rcpr_D_squared, cos_B, cos_phi, cos_phi0;
    long double cos_phi_phi0, psi;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0L / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B        = rssringoccs_LDouble_Cos(B);
    cos_phi      = rssringoccs_LDouble_Cos(phi);
    cos_phi0     = rssringoccs_LDouble_Cos(phi0);
    cos_phi_phi0 = rssringoccs_LDouble_Cos(phi-phi0);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos_B * rcpr_D) * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r - 2.0L*r*r0*cos_phi_phi0) * rcpr_D_squared;

    /* Sign of xi swapped from MTR86.                                         */
    psi = kD * (rssringoccs_LDouble_Sqrt(1.0L + eta - 2.0L*xi) + xi - 1.0L);

    return psi;
}
