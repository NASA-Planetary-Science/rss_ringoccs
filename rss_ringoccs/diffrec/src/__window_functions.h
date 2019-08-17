#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

/*  Various coefficients and constants defined here.    */
#include "__math_constants.h"

/*  Kaiser-Bessel Functions found here.                 */
#include "__kaiser_bessel.h"
#include "__modified_kaiser_bessel.h"

float Rect_Window_Float(float x, float W)
{
    if (fabsf(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

double Rect_Window_Double(double x, double W)
{
    if (fabs(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

long double Rect_Window_Long_Double(long double x, long double W)
{
    if (fabsl(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

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

#endif