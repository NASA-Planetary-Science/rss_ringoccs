#ifndef RSS_RINGOCCS_SINC_H
#define RSS_RINGOCCS_SINC_H

#include <math.h>

float Sinc_Float(float x){
    if (x == 0){
        return 1.0;
    }
    else{
        return sinf(x)/x;
    }
}

double Sinc_Double(double x){
    if (x == 0){
        return 1.0;
    }
    else{
        return sin(x)/x;
    }
}

long double Sinc_Long_Double(long double x){
    if (x == 0){
        return 1.0;
    }
    else{
        return sinl(x)/x;
    }
}

#endif