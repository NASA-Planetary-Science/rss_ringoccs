#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double Newton_Raphson_Fresnel_Psi_D(double k, double r, double r0,
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
    dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi(k, r, r0, phi, phi0, B, D);

    while(fabs(dphi) > EPS)
    {
        dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi(k, r, r0, phi, phi0, B, D);
        phi  -= dphi/rssringoccs_Double_Fresnel_d2Psi_dPhi2(k, r, r0, phi,
                                                            phi0, B, D);
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
