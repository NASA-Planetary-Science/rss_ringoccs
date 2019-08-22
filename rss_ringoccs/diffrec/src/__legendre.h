#include <stdio.h>

void Legendre_Coeffs(double *poly, double x, int order){
    int i;
    poly[0] = 1.0;
    poly[1] = x;
    for (i=2; i<order; ++i){
        poly[i] = ((2.0*i-1.0)*x*poly[i-1]-(i-1.0)*poly[i-2])/i;
    }
}

void Alt_Legendre_Coeffs(double *poly, double *legendre_p, int order){
    int i;
    for (i=0; i<order; ++i){
        poly[i] = (legendre_p[i]-legendre_p[1]*legendre_p[i+1])/(i+2.0);
    }
}