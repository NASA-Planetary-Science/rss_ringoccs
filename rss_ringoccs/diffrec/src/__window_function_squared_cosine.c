#include "__window_functions.h"

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
double Coss_Window_Char(char x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_UChar(unsigned char x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_Short(short x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_UShort(unsigned short x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_Int(int x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_UInt(unsigned int x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_Long(long x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_ULong(unsigned long x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_Long_Long(long long x, double W)
{
    return Coss_Window_Double((double)x, W);
}

double Coss_Window_ULong_Long(unsigned long long x, double W)
{
    return Coss_Window_Double((double)x, W);
}