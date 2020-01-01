#include "__math_functions.h"

float Sinc_Float(float x){
    if (x == 0) return 1.0;
    else return sinf(x)/x;
}

double Sinc_Double(double x){
    if (x == 0) return 1.0;
    else return sin(x)/x;
}

long double Sinc_Long_Double(long double x){
    if (x == 0) return 1.0;
    else return sinl(x)/x;
}

/*  For all integer types, convert to double and compute.                     */
double Sinc_Char(char x)                     {return Sinc_Double((double)x);}
double Sinc_UChar(unsigned char x)           {return Sinc_Double((double)x);}
double Sinc_Short(short x)                   {return Sinc_Double((double)x);}
double Sinc_UShort(unsigned short x)         {return Sinc_Double((double)x);}
double Sinc_Int(int x)                       {return Sinc_Double((double)x);}
double Sinc_UInt(unsigned int x)             {return Sinc_Double((double)x);}
double Sinc_Long(long x)                     {return Sinc_Double((double)x);}
double Sinc_ULong(unsigned long x)           {return Sinc_Double((double)x);}
double Sinc_Long_Long(long long x)           {return Sinc_Double((double)x);}
double Sinc_ULong_Long(unsigned long long x) {return Sinc_Double((double)x);}