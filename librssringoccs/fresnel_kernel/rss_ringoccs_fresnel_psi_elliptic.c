#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float
rssringoccs_Float_Fresnel_dPsi_dPhi_Ellipse(float kD,  float r,    float r0,
                                            float phi, float phi0, float B,
                                            float D,   float ecc,  float peri)
{
    float xi, eta, psi0, psi_d1, dxi_phi, deta_phi, dxi_rho, deta_rho;
    float cos_B, cos_phi, sin_phi, cos_phi0, sin_phi0;
    float rcpr_D, rcpr_D_squared, cos_phi_phi0, sin_phi_phi0;
    float xi_factor, eta_factor, dpsi_rho, dpsi_phi, drho_phi, semi_major;
    float cos_phi_peri, sin_phi_peri, ecc_factor, ecc_cos_factor;

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

    cos_phi_peri = rssringoccs_Float_Cos(phi - peri);
    sin_phi_peri = rssringoccs_Float_Sin(phi - peri);

    /*  These terms appear do to the eccentricity.                            */
    ecc_factor     = 1.0F - ecc*ecc;
    ecc_cos_factor = 1.0F + ecc * cos_phi_peri;
    semi_major     = r * ecc_cos_factor / ecc_factor;

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  This term appears in dxi and dxi2 and xi.                             */
    xi_factor = cos_B * rcpr_D;

    /*  And this term appears in eta, deta, and deta2.                        */
    eta_factor = 2.0F * r * r0 * rcpr_D_squared;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = xi_factor * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - eta_factor*cos_phi_phi0;
    psi0 = rssringoccs_Float_Sqrt(1.0F + eta - 2.0F*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -xi_factor * (r * sin_phi);
    deta_phi = eta_factor * sin_phi_phi0;
    dxi_rho  = xi_factor * cos_phi;
    deta_rho = 2.0F * (r - r0 * cos_phi_phi0) * rcpr_D_squared;

    drho_phi = semi_major*ecc*ecc_factor*sin_phi_peri /
               (ecc_cos_factor*ecc_cos_factor);

    dpsi_rho = (deta_rho - 2.0F * dxi_rho) * (0.5F / psi0) + dxi_rho;
    dpsi_phi = (deta_phi - 2.0F * dxi_phi) * (0.5F / psi0) + dxi_phi;

    /*  Compute the partial derivative.                                       */
    psi_d1  = kD*(dpsi_phi + dpsi_rho*drho_phi);
    return psi_d1;
}

double
rssringoccs_Double_Fresnel_dPsi_dPhi_Ellipse(double kD, double r, double r0,
                                             double phi, double phi0, double B,
                                             double D, double ecc, double peri)
{
    double xi, eta, psi0, psi_d1, dxi_phi, deta_phi, dxi_rho, deta_rho;
    double cos_B, cos_phi, sin_phi, cos_phi0, sin_phi0;
    double rcpr_D, rcpr_D_squared, cos_phi_phi0, sin_phi_phi0;
    double xi_factor, eta_factor, dpsi_rho, dpsi_phi, drho_phi, semi_major;
    double cos_phi_peri, sin_phi_peri, ecc_factor, ecc_cos_factor;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0 / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B        = rssringoccs_Double_Cos(B);
    cos_phi      = rssringoccs_Double_Cos(phi);
    sin_phi      = rssringoccs_Double_Sin(phi);
    cos_phi0     = rssringoccs_Double_Cos(phi0);
    sin_phi0     = rssringoccs_Double_Sin(phi0);
    cos_phi_peri = rssringoccs_Double_Cos(phi - peri);
    sin_phi_peri = rssringoccs_Double_Sin(phi - peri);

    /*  These terms appear do to the eccentricity.                            */
    ecc_factor     = 1.0 - ecc*ecc;
    ecc_cos_factor = 1.0 + ecc * cos_phi_peri;
    semi_major     = r * ecc_cos_factor / ecc_factor;

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
    psi0 = rssringoccs_Double_Sqrt(1.0 + eta - 2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -xi_factor * (r * sin_phi);
    deta_phi = eta_factor * sin_phi_phi0;
    dxi_rho  = xi_factor * cos_phi;
    deta_rho = 2.0 * (r - r0 * cos_phi_phi0) * rcpr_D_squared;

    drho_phi = semi_major*ecc*ecc_factor*sin_phi_peri /
               (ecc_cos_factor*ecc_cos_factor);

    dpsi_rho = (deta_rho - 2.0 * dxi_rho) * (0.5 / psi0) + dxi_rho;
    dpsi_phi = (deta_phi - 2.0 * dxi_phi) * (0.5 / psi0) + dxi_phi;

    /*  Compute the partial derivative.                                       */
    psi_d1  = kD*(dpsi_phi + dpsi_rho*drho_phi);
    return psi_d1;
}

long double
rssringoccs_LDouble_Fresnel_dPsi_dPhi_Ellipse(long double kD,   long double r,
                                              long double r0,   long double phi,
                                              long double phi0, long double B,
                                              long double D,    long double ecc,
                                              long double peri)
{
    long double xi, eta, psi0, psi_d1, dxi_phi, deta_phi, dxi_rho, deta_rho;
    long double cos_B, cos_phi, sin_phi, cos_phi0, sin_phi0;
    long double rcpr_D, rcpr_D_squared, cos_phi_phi0, sin_phi_phi0;
    long double xi_factor, eta_factor, dpsi_rho, dpsi_phi, drho_phi, semi_major;
    long double cos_phi_peri, sin_phi_peri, ecc_factor, ecc_cos_factor;

    /*  Compute 1/D and it's square to save the number of divisions we need   *
     *  to compute. Multiplication is usually ~10 times faster.               */
    rcpr_D = 1.0L / D;
    rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Precompute cosines and sines to save on computations.                 */
    cos_B        = rssringoccs_LDouble_Cos(B);
    cos_phi      = rssringoccs_LDouble_Cos(phi);
    sin_phi      = rssringoccs_LDouble_Sin(phi);
    cos_phi0     = rssringoccs_LDouble_Cos(phi0);
    sin_phi0     = rssringoccs_LDouble_Sin(phi0);
    cos_phi_peri = rssringoccs_LDouble_Cos(phi - peri);
    sin_phi_peri = rssringoccs_LDouble_Sin(phi - peri);

    /*  These terms appear do to the eccentricity.                            */
    ecc_factor     = 1.0L - ecc*ecc;
    ecc_cos_factor = 1.0L + ecc * cos_phi_peri;
    semi_major     = r * ecc_cos_factor / ecc_factor;

    /*  Since we've computed cos and sin of phi and phi0, cos and sin of      *
     *  phi-phi0 can be computed without the need to call cos and sin again.  */
    cos_phi_phi0 = sin_phi*sin_phi0 + cos_phi*cos_phi0;
    sin_phi_phi0 = sin_phi*cos_phi0 - cos_phi*sin_phi0;

    /*  This term appears in dxi and dxi2 and xi.                             */
    xi_factor = cos_B * rcpr_D;

    /*  And this term appears in eta, deta, and deta2.                        */
    eta_factor = 2.0L * r * r0 * rcpr_D_squared;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = xi_factor * (r * cos_phi - r0 * cos_phi0);
    eta  = (r0*r0 + r*r)*rcpr_D_squared - eta_factor*cos_phi_phi0;
    psi0 = rssringoccs_LDouble_Sqrt(1.0 + eta - 2.0L*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -xi_factor * (r * sin_phi);
    deta_phi = eta_factor * sin_phi_phi0;
    dxi_rho  = xi_factor * cos_phi;
    deta_rho = 2.0L * (r - r0 * cos_phi_phi0) * rcpr_D_squared;

    drho_phi = semi_major*ecc*ecc_factor*sin_phi_peri /
               (ecc_cos_factor*ecc_cos_factor);

    dpsi_rho = (deta_rho - 2.0L * dxi_rho) * (0.5L / psi0) + dxi_rho;
    dpsi_phi = (deta_phi - 2.0L * dxi_phi) * (0.5L / psi0) + dxi_phi;

    /*  Compute the partial derivative.                                       */
    psi_d1  = kD*(dpsi_phi + dpsi_rho*drho_phi);
    return psi_d1;
}
