#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double
rssringoccs_Double_Newton_Raphson_Fresnel_Ellipse(double kD, double r,
                                                  double r0, double phi,
                                                  double phi0, double B,
                                                  double D, double ecc,
                                                  double peri, double EPS,
                                                  unsigned char toler)
{
    double dphi, err, semi_major, ecc_factor, ecc_cos_factor, cos_phi_peri;
    unsigned char n = 0;
    dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi_Ellipse(kD, r, r0, phi, phi0,
                                                         B, D, ecc, peri);

    dphi /= rssringoccs_Double_Fresnel_d2Psi_dPhi2(kD, r, r0, phi,
                                                   phi0, B, D);

    cos_phi_peri = rssringoccs_Double_Cos(phi - peri);
    ecc_factor     = 1.0 - ecc*ecc;
    ecc_cos_factor = 1.0 + ecc * cos_phi_peri;
    semi_major     = r * ecc_cos_factor / ecc_factor;
    phi  -= dphi;
    err = rssringoccs_Double_Abs(dphi);
    while(err > EPS)
    {
        cos_phi_peri = rssringoccs_Double_Cos(phi - peri);
        ecc_cos_factor = 1.0 + ecc * cos_phi_peri;
        r = semi_major * ecc_factor / (1.0 + ecc*cos_phi_peri);

        dphi  = rssringoccs_Double_Fresnel_dPsi_dPhi_Ellipse(kD, r, r0, phi,
                                                             phi0, B, D,
                                                             ecc, peri);

        dphi /= rssringoccs_Double_Fresnel_d2Psi_dPhi2(kD, r, r0, phi,
                                                       phi0, B, D);
        phi  -= dphi;
        ++n;
        if (n > toler)
            break;

        err = rssringoccs_Double_Abs(dphi);
    }
    return phi;
}