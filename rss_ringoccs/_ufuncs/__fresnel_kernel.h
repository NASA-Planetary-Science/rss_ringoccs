
/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_H

#include <stdio.h>
#include <math.h>

double Newton_Raphson(double x, double (*f)(double),
                      double (*f_prime)(double), double EPS, long toler)
{
    double dx;
    long i = 0;
    dx = (*f)(x)/(*f_prime)(x);
    x -= dx;
    while(fabs(dx) > EPS){
        dx = (*f)(x)/(*f_prime)(x);
        x -= dx;
        ++i;
        if (i > toler){break;}
    }
    return x;
}

double fresnel_psi(double kD, double r, double r0, double phi,
                   double phi0, double B, double D)
{
    double xi, eta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrt(1.0+eta-2.0*xi) + xi - 1.0);
}

double fresnel_dpsi_dphi(double kD, double r, double r0, double phi, 
                         double phi0, double B, double D)
{
    double xi, eta, psi0, dxi, deta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives (Use calculus before hand).                       */
    dxi  = -(cos(B) / D) * (r * sin(phi));
    deta = 2.0 * r * r0 * sin(phi - phi0) / (D*D);

    return kD * ((0.5/psi0)*(deta-2.0*dxi) + dxi);
}

double fresnel_dpsi_dphi_ellipse(double kD, double r, double r0,
                                 double phi, double phi0, double B,
                                 double D, double ecc, double peri)
{
    double xi, eta, psi0, psi_d1, dxi_phi, deta_phi, dxi_rho, deta_rho;
    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);
    psi0 = sqrt(1.0+eta-2.0*xi);

    /*  Compute derivatives.                                                  */
    dxi_phi  = -(cos(B) / D) * (r * sin(phi));
    deta_phi = 2.0 * r * r0 * sin(phi-phi0) / (D*D);
    dxi_rho  = (cos(B) / D) * cos(phi);
    deta_rho = 2.0 * (r - r0 * cos(phi - phi0)) / (D*D);

    /*  Compute the partial derivative.                                       */
    psi_d1  = (deta_rho-2.0*dxi_rho)*(0.5/psi0) + dxi_rho;
    psi_d1 *= r * ecc * sin(phi - peri) / (1 + ecc * cos(phi - peri));
    psi_d1 += (deta_phi - 2.0 * dxi_phi) * (0.5 / psi0) + dxi_phi;

    return kD * psi_d1;
}

double fresnel_d2psi_dphi2(double kD, double r, double r0, double phi,
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

double Newton_Raphson_fresnel_psi(double kD, double r, double r0,
                                  double phi, double phi0, double B,
                                  double D, double EPS, long toler)
{
    double dphi;
    long i = 0;
    dphi  = fresnel_psi(kD, r, r0, phi, phi0, B, D);
    while(fabs(dphi) > EPS){
        dphi  = fresnel_dpsi_dphi(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi/fresnel_d2psi_dphi2(kD, r, r0, phi, phi0, B, D);
        ++i;
        if (i > toler){break;}
    }
    return phi;
}

double Newton_Raphson_fresnel_ellipse(double kD, double r, double r0,
                                      double phi, double phi0, double B,
                                      double D, double ecc, double peri,
                                      double EPS, long toler)
{
    double dphi;
    long i = 0;
    dphi  = fresnel_dpsi_dphi_ellipse(kD, r, r0, phi, phi0, B, D, ecc, peri);
    dphi /= fresnel_d2psi_dphi2(kD, r, r0, phi, phi0, B, D);
    phi  -= dphi;
    while(fabs(dphi) > EPS){
        dphi  = fresnel_dpsi_dphi_ellipse(kD, r, r0, phi, phi0, B, D, ecc, peri);
        dphi /= fresnel_d2psi_dphi2(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi;
        ++i;
        if (i > toler){break;}
    }
    return phi;
}

#endif