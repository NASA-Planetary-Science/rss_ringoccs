/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_H

/*  Various mathematical functions defined here.                              */
#include <math.h>

/*---------------------------Newton-Raphson Function--------------------------*/

float Newton_Raphson_Float(float x, float (*f)(float),
                           float (*f_prime)(float), float EPS, long toler){
    float dx;
    long i = 0;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    dx = (*f)(x)/(*f_prime)(x);
    x -= dx;

    /*  Continuing this computation until the error is below the threshold.   */
    while(fabsf(dx) > EPS){
        dx = (*f)(x)/(*f_prime)(x);
        x -= dx;
        ++i;

        /*  Break if too many iterations have been run.                       */
        if (i > toler){
            break;
        }
    }
    return x;
}

double Newton_Raphson_Double(double x, double (*f)(double),
                             double (*f_prime)(double), double EPS, long toler){
    double dx;
    long i = 0;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    dx = (*f)(x)/(*f_prime)(x);
    x -= dx;

    /*  Continuing this computation until the error is below the threshold.   */
    while(fabs(dx) > EPS){
        dx = (*f)(x)/(*f_prime)(x);
        x -= dx;
        ++i;

        /*  Break if too many iterations have been run.                       */
        if (i > toler){
            break;
        }
    }
    return x;
}

long double Newton_Raphson_Long_Double(long double x,
                                       long double (*f)(long double),
                                       long double (*f_prime)(long double),
                                       long double EPS, long toler){
    long double dx;
    long i = 0;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    dx = (*f)(x)/(*f_prime)(x);
    x -= dx;

    /*  Continuing this computation until the error is below the threshold.   */
    while(fabsl(dx) > EPS){
        dx = (*f)(x)/(*f_prime)(x);
        x -= dx;
        ++i;

        /*  Break if too many iterations have been run.                       */
        if (i > toler){
            break;
        }
    }
    return x;
}

/*----------------------------The Fresnel Kernel------------------------------*/

float Fresnel_Psi_Float(float kD, float r, float r0, float phi,
                        float phi0, float B, float D){
    float xi, eta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cosf(B)/D) * (r * cosf(phi) - r0 * cosf(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cosf(phi-phi0)) / (D*D);

    /* Sign of xi swapped from MTR86.                                         */
    return kD * (sqrtf(1.0+eta-2.0*xi) + xi - 1.0);
}

double Fresnel_Psi_Double(double kD, double r, double r0, double phi,
                          double phi0, double B, double D){
    double xi, eta;

    /*  Compute xi variable (MTR86 Equation 4b) and eta (Equation 4c).        */
    xi   = (cos(B)/D) * (r * cos(phi) - r0 * cos(phi0));
    eta  = (r0*r0 + r*r - 2.0*r*r0*cos(phi-phi0)) / (D*D);

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

/*----------------The First Derivative of the Fresnel Kernel------------------*/

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

/*---------------The Second Derivative of the Fresnel Kernel------------------*/

double Fresnel_dPsi_dPhi_Func(double kD, double r, double r0, double phi,
                              double phi0, double B, double D){
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

double Fresnel_dPsi_dPhi_Ellipse_Func(double kD, double r, double r0,
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

double Fresnel_d2Psi_dPhi2_Func(double kD, double r, double r0, double phi,
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

double Newton_Raphson_Fresnel_Psi(double kD, double r, double r0,
                                  double phi, double phi0, double B,
                                  double D, double EPS, long toler)
{
    double dphi;
    long i = 0;
    dphi  = Fresnel_dPsi_dPhi_Double(kD, r, r0, phi, phi0, B, D);
    while(fabs(dphi) > EPS){
        dphi  = Fresnel_dPsi_dPhi_Double(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi/Fresnel_d2Psi_dPhi2_Func(kD, r, r0, phi, phi0, B, D);
        ++i;
        if (i > toler){
            break;
        }
    }
    return phi;
}

double Newton_Raphson_Fresnel_Ellipse(double kD, double r, double r0,
                                      double phi, double phi0, double B,
                                      double D, double ecc, double peri,
                                      double EPS, long toler)
{
    double dphi;
    long i = 0;
    dphi  = Fresnel_dPsi_dPhi_Ellipse_Func(kD, r, r0, phi, phi0,
                                           B, D, ecc, peri);
    dphi /= Fresnel_d2Psi_dPhi2_Func(kD, r, r0, phi, phi0, B, D);
    phi  -= dphi;
    while(fabs(dphi) > EPS){
        dphi  = Fresnel_dPsi_dPhi_Ellipse_Func(kD, r, r0, phi, phi0,
                                               B, D, ecc, peri);
        dphi /= Fresnel_d2Psi_dPhi2_Func(kD, r, r0, phi, phi0, B, D);
        phi  -= dphi;
        ++i;
        if (i > toler){break;}
    }
    return phi;
}

#endif
