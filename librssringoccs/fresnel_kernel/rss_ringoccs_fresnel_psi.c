#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

float Fresnel_Psi_Float(float kD, float r, float r0, float phi,
                        float phi0, float B, float D)
{
    float xi, eta;
    float rcpr_D = 1.0 / D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosf(B) * rcpr_D) * (r * cosf(phi) - r0 * cosf(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosf(phi-phi0)) * rcpr_D * rcpr_D;

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrtf(1.0+eta-2.0*xi) + xi - 1.0);
}

double Fresnel_Psi_Double(double kD, double r, double r0, double phi,
                          double phi0, double B, double D){
    double xi, eta;
    double rcpr_D = 1.0 / D;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B) * rcpr_D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) * rcpr_D * rcpr_D;

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrt(1.0+eta-2.0*xi) + xi - 1.0);
}

double Fresnel_Psi_Long_Double(long double kD, long double r, long double r0,
                               long double phi, long double phi0, long double B,
                               long double D){
    long double xi, eta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosl(B)/D) * (r * cosl(phi) - r0 * cosl(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrtl(1.0+eta-2.0*xi) + xi - 1.0);
}
