#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float Fresnel_dPsi_dPhi_Ellipse_Float(float kD,  float r,    float r0,
                                      float phi, float phi0, float B,
                                      float D,   float ecc,  float peri)
{
    float xi, eta, psi0, psi_d1;
    float dxi_phi, deta_phi, dxi_rho, deta_rho;
    float rcpr_D = 1.0/D;
    float rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosf(B) * rcpr_D) * (r * cosf(phi) - r0 * cosf(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosf(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtf(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -(cosf(B) * rcpr_D) * (r * sinf(phi));
    deta_phi = 2.0 * r * r0 * sinf(phi-phi0) * rcpr_D_squared;
    dxi_rho  = (cosf(B) * rcpr_D) * cosf(phi);
    deta_rho = 2.0 * (r - r0 * cosf(phi - phi0)) * rcpr_D_squared;

    /*  Compute the partial derivative.                                       */
    psi_d1  = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho;
    psi_d1 *= r * ecc * sinf(phi - peri) / (1 + ecc * cosf(phi - peri));
    psi_d1 += (deta_phi - 2.0 * dxi_phi) * (0.5 / psi0) + dxi_phi;

    return kD * psi_d1;
}

double Fresnel_dPsi_dPhi_Ellipse_Double(double kD,  double r,    double r0,
                                        double phi, double phi0, double B,
                                        double D,   double ecc,  double peri)
{
    double xi, eta, psi0, psi_d1;
    double dxi_phi, deta_phi, dxi_rho, deta_rho;
    double rcpr_D = 1.0/D;
    double rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B) * rcpr_D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -(cos(B) * rcpr_D) * (r * sin(phi));
    deta_phi = 2.0 * r * r0 * sin(phi-phi0) * rcpr_D_squared;
    dxi_rho  = (cos(B) * rcpr_D) * cos(phi);
    deta_rho = 2.0 * (r - r0 * cos(phi - phi0)) * rcpr_D_squared;

    /*  Compute the partial derivative.                                       */
    psi_d1  = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho;
    psi_d1 *= r * ecc * sin(phi - peri) / (1 + ecc * cos(phi - peri));
    psi_d1 += (deta_phi - 2.0 * dxi_phi) * (0.5 / psi0) + dxi_phi;

    return kD * psi_d1;
}

long double Fresnel_dPsi_dPhi_Ellipse_Long_Double(
    long double kD,   long double r, long double r0, long double phi,
    long double phi0, long double B, long double D,  long double ecc,
    long double peri
)
{
    long double xi, eta, psi0, psi_d1;
    long double dxi_phi, deta_phi, dxi_rho, deta_rho;
    long double rcpr_D = 1.0/D;
    long double rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosl(B) * rcpr_D) * (r * cosl(phi) - r0 * cosl(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosl(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtl(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -(cosl(B) * rcpr_D) * (r * sinl(phi));
    deta_phi = 2.0 * r * r0 * sinl(phi-phi0) * rcpr_D_squared;
    dxi_rho  = (cosl(B) * rcpr_D) * cosl(phi);
    deta_rho = 2.0 * (r - r0 * cosl(phi - phi0)) * rcpr_D_squared;

    /*  Compute the partial derivative.                                       */
    psi_d1  = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho;
    psi_d1 *= r * ecc * sinl(phi - peri) / (1 + ecc * cosl(phi - peri));
    psi_d1 += (deta_phi - 2.0 * dxi_phi) * (0.5 / psi0) + dxi_phi;

    return kD * psi_d1;
}
