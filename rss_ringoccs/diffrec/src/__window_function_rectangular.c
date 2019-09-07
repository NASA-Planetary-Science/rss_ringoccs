#include "__window_functions.h"

float Rect_Window_Float(float x, float W){
    if (fabsf(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

double Rect_Window_Double(double x, double W){
    if (fabs(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

long double Rect_Window_Long_Double(long double x, long double W){
    if (fabsl(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}