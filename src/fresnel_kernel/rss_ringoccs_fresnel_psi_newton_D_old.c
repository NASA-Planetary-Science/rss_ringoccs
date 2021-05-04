#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

static double __dpsi(double kD, double r, double r0, double phi,
                     double phi0, double B, double D)
{
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

static double __d2psi(double kD, double r, double r0, double phi,
                      double phi0, double B, double D)
{
    double xi, eta, psi0, dxi, dxi2, deta, deta2, psi_d2;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi   = -(cos(B)/D) * (r * sin(phi));
    dxi2  = -(cos(B)/D) * (r * cos(phi));
    deta  = 2.0 * r * r0 * sin(phi - phi0) / (D*D);
    deta2 = 2.0 * r * r0 * cos(phi - phi0) / (D*D);

    /*  Compute the second partial derivative.                                */
    psi_d2  = (-0.25/(psi0*psi0*psi0)) * (deta - 2.0*dxi) * (deta - 2.0*dxi);
    psi_d2 += (0.5/psi0) * (deta2 - 2.0*dxi2) + dxi2;

    return kD*psi_d2;
}

double Newton_Raphson_Fresnel_Psi_D_Old(double kD, double r, double r0,
                                        double phi, double phi0, double B,
                                        double EPS, long toler, double rx,
                                        double ry, double rz)
{
    double dphi;
    double D;
    double x, y, dx, dy;
    long i = 0;
    x = r0 * cos(phi);
    y = r0 * sin(phi);
    dx = x-rx;
    dy = y-ry;
    D = sqrt(dx*dx + dy*dy + rz*rz);
    dphi  = __dpsi(kD, r, r0, phi, phi0, B, D);

    while(fabs(dphi) > EPS){
        dphi  = __dpsi(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi/__d2psi(kD, r, r0, phi, phi0, B, D);
        ++i;
        if (i > toler)
            break;

        x = r0 * cos(phi);
        y = r0 * sin(phi);
        dx = x-rx;
        dy = y-ry;
        D = sqrt(dx*dx + dy*dy + rz*rz);
    }
    return phi;
}
