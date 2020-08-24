#include "special_functions.h"
#include "_math_constants.h"
#include <math.h>

float Coss_Window_Float(float x, float W)
{
    if (fabsf(x) <= W/2.0){
        x *= ONE_PI/W;
        x = cos(x);
        x *=x;
        return x;
    }
    else {
        return 0.0;
    }
}

double Coss_Window_Double(double x, double W)
{
    if (fabs(x) <= W/2.0){
        x *= ONE_PI/W;
        x = cos(x);
        x *=x;
        return x;
    }
    else {
        return 0.0;
    }
}

long double Coss_Window_Long_Double(long double x, long double W)
{
    if (fabsl(x) <= W/2.0){
        x *= ONE_PI/W;
        x = cos(x);
        x *=x;
        return x;
    }
    else {
        return 0.0;
    }
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputTwoVarForFloatOutput(Coss_Window);