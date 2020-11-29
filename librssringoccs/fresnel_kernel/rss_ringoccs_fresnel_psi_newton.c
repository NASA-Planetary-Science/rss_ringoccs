#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

double Newton_Raphson_Fresnel_Psi(double kD, double r, double r0,
                                  double phi, double phi0, double B,
                                  double D, double EPS, long toler){
    double dphi;
    long i = 0;
    dphi  = Fresnel_dPsi_dPhi_Double(kD, r, r0, phi, phi0, B, D);
    while(fabs(dphi) > EPS){
        dphi  = Fresnel_dPsi_dPhi_Double(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi/Fresnel_d2Psi_dPhi2_Double(kD, r, r0, phi, phi0, B, D);
        ++i;
        if (i > toler)
            break;
    }
    return phi;
}

double Newton_Raphson_Fresnel_Ellipse(double kD, double r, double r0,
                                      double phi, double phi0, double B,
                                      double D, double ecc, double peri,
                                      double EPS, unsigned char toler){
    double dphi;
    unsigned char i = 0;
    dphi  = Fresnel_dPsi_dPhi_Ellipse_Double(kD, r, r0, phi, phi0,
                                             B, D, ecc, peri);
    while(fabs(dphi) > EPS){
        dphi  = Fresnel_dPsi_dPhi_Ellipse_Double(kD, r, r0, phi, phi0,
                                                 B, D, ecc, peri);
        dphi /= Fresnel_d2Psi_dPhi2_Double(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi;
        ++i;
        if (i > toler){
            break;
        }
    }
    return phi;
}
