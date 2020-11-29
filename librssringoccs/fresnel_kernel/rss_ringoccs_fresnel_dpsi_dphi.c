#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float Fresnel_dPsi_dPhi_Float(float kD, float r, float r0, float phi,
                              float phi0, float B, float D){
    float xi, eta, psi0, dxi, deta;

    float cosB_over_D = cosf(B)/D;
    float rcpr_D_squared = 1.0 / (D*D);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cosB_over_D * (r * cosf(phi) - r0 * cosf(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosf(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtf(1.0+eta-2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cosB_over_D * (r * sinf(phi));
    deta = 2.0 * r * r0 * sinf(phi - phi0) * rcpr_D_squared;

    return kD * ((0.5/psi0)*(deta-2.0*dxi) + dxi);
}

double Fresnel_dPsi_dPhi_Double(double kD, double r, double r0, double phi,
                                double phi0, double B, double D){
    double xi, eta, psi0, dxi, deta;

    double cosB_over_D = cos(B)/D;
    double rcpr_D_squared = 1.0 / (D*D);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cosB_over_D * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cosB_over_D * (r * sin(phi));
    deta = 2.0 * r * r0 * sin(phi - phi0) * rcpr_D_squared;

    return kD * ((0.5/psi0)*(deta-2.0*dxi) + dxi);
}

double Fresnel_dPsi_dPhi_Long_Double(long double kD, long double r,
                                     long double r0, long double phi,
                                     long double phi0, long double B,
                                     long double D){
    long double xi, eta, psi0, dxi, deta;

    long double cosB_over_D = cosl(B)/D;
    long double rcpr_D_squared = 1.0 / (D*D);

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = cosB_over_D * (r * cosl(phi) - r0 * cosl(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosl(phi-phi0)) * rcpr_D_squared;
    psi0 = sqrtl(1.0+eta-2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -cosB_over_D * (r * sinl(phi));
    deta = 2.0 * r * r0 * sinl(phi - phi0) * rcpr_D_squared;

    return kD * ((0.5/psi0)*(deta-2.0*dxi) + dxi);
}
