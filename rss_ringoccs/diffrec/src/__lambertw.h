/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_LAMBERTW_H
#define RSS_RINGOCCS_LAMBERTW_H

#include "__math_constants.h"
#include <math.h>

float LambertW_Float(float x){
    float EPS = 1.0e-4;
    float x0, dx, exp_x0;
    if (x > 2.0){
        x0 = logf(x/logf(x));
    }
    else if (x > -RCPR_EULER_E){
        x0 = x;
    }
    else if (x == -RCPR_EULER_E){
        return -1.0;
    }
    else {
        return nanf("0");
    }
    exp_x0 = expf(x0);
    dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
    x0 = x0 - dx;
    while (fabsf(dx) > EPS){
        exp_x0 = expf(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;
    }
    return x0;
}

double LambertW_Double(double x){
    double EPS = 1.0e-8;
    double x0, dx, exp_x0;
    if (x > 2.0){
        x0 = log(x/log(x));
    }
    else if (x > -RCPR_EULER_E){
        x0 = x;
    }
    else if (x == -RCPR_EULER_E){
        return -1.0;
    }
    else {
        return nan("0");
    }
    exp_x0 = exp(x0);
    dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
    x0 = x0 - dx;
    while (fabs(dx) > EPS){
        exp_x0 = exp(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;
    }
    return x0;
}

long double LambertW_Long_Double(long double x){
    long double EPS = 1.0e-16;
    long double x0, dx, exp_x0;
    if (x > 2.0){
        x0 = logl(x/logl(x));
    }
    else if (x > -RCPR_EULER_E){
        x0 = x;
    }
    else if (x == -RCPR_EULER_E){
        return -1.0;
    }
    else {
        return nanl("0");
    }
    exp_x0 = expl(x0);
    dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
    x0 = x0 - dx;
    while (fabsl(dx) > EPS){
        exp_x0 = expl(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0)-(x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;
    }
    return x0;
}

#endif