#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

/*  Various coefficients and constants defined here.    */
#include "__math_constants.h"

/*  Kaiser-Bessel Functions found here.                 */
#include "__kaiser_bessel.h"
#include "__modified_kaiser_bessel.h"

double Rect_Window_Func(double x, double W)
{
    if (fabs(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

double Coss_Window_Func(double x, double W)
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

#endif