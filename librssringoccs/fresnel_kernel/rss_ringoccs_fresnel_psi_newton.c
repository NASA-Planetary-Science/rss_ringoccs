#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double Newton_Raphson_Fresnel_Psi(double kD, double r, double r0,
                                  double phi, double phi0, double B,
                                  double D, double EPS, unsigned char toler)
{
    double dphi, err;
    unsigned char n = 0;
    dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi(kD, r, r0, phi, phi0, B, D);
    err = rssringoccs_Double_Abs(dphi);
    while(err > EPS)
    {
        dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi(kD, r, r0, phi,
                                                     phi0, B, D);
        phi  -= dphi/rssringoccs_Double_Fresnel_d2Psi_dPhi2(kD, r, r0, phi,
                                                            phi0, B, D);
        ++n;
        if (n > toler)
            break;

        err = rssringoccs_Double_Abs(dphi);
    }
    return phi;
}
