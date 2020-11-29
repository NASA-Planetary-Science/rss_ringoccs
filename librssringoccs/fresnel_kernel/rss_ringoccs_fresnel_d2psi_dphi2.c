#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float Fresnel_d2Psi_dPhi2_Float(float kD,   float r, float r0, float phi,
                                float phi0, float B, float D)
{
    float xi, eta, psi0, dxi, dxi2, deta, deta2, psi_d2;
    float rcpr_D = 1.0 / D;
    float rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosf(B) *  rcpr_D) * (r * cosf(phi) - r0 * cosf(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosf(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtf(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi   = -(cosf(B) * rcpr_D) * (r * sinf(phi));
    dxi2  = -(cosf(B) * rcpr_D) * (r * cosf(phi));
    deta  = 2.0 * r * r0 * sinf(phi - phi0) * rcpr_D_squared;
    deta2 = 2.0 * r * r0 * cosf(phi - phi0) * rcpr_D_squared;

    /*  Compute the second partial derivative.                                */
    psi_d2  = (-0.25/(psi0*psi0*psi0)) * (deta - 2.0*dxi) * (deta - 2.0*dxi);
    psi_d2 += (0.5/psi0) * (deta2 - 2.0*dxi2) + dxi2;

    return kD*psi_d2;
}

double Fresnel_d2Psi_dPhi2_Double(double kD,   double r, double r0, double phi,
                                  double phi0, double B, double D)
{
    double xi, eta, psi0, dxi, dxi2, deta, deta2, psi_d2;
    double rcpr_D = 1.0 / D;
    double rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B) *  rcpr_D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi   = -(cos(B) * rcpr_D) * (r * sin(phi));
    dxi2  = -(cos(B) * rcpr_D) * (r * cos(phi));
    deta  = 2.0 * r * r0 * sin(phi - phi0) * rcpr_D_squared;
    deta2 = 2.0 * r * r0 * cos(phi - phi0) * rcpr_D_squared;

    /*  Compute the second partial derivative.                                */
    psi_d2  = (-0.25/(psi0*psi0*psi0)) * (deta - 2.0*dxi) * (deta - 2.0*dxi);
    psi_d2 += (0.5/psi0) * (deta2 - 2.0*dxi2) + dxi2;

    return kD*psi_d2;
}

long double Fresnel_d2Psi_dPhi2_Long_Double(long double kD,   long double r,
                                            long double r0,   long double phi,
                                            long double phi0, long double B,
                                            long double D)
{
    long double xi, eta, psi0, dxi, dxi2, deta, deta2, psi_d2;
    long double rcpr_D = 1.0 / D;
    long double rcpr_D_squared = rcpr_D * rcpr_D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosl(B) *  rcpr_D) * (r * cosl(phi) - r0 * cosl(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosl(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtl(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi   = -(cosl(B) * rcpr_D) * (r * sinl(phi));
    dxi2  = -(cosl(B) * rcpr_D) * (r * cosl(phi));
    deta  = 2.0 * r * r0 * sinl(phi - phi0) * rcpr_D_squared;
    deta2 = 2.0 * r * r0 * cosl(phi - phi0) * rcpr_D_squared;

    /*  Compute the second partial derivative.                                */
    psi_d2  = (-0.25/(psi0*psi0*psi0)) * (deta - 2.0*dxi) * (deta - 2.0*dxi);
    psi_d2 += (0.5/psi0) * (deta2 - 2.0*dxi2) + dxi2;

    return kD*psi_d2;
}
