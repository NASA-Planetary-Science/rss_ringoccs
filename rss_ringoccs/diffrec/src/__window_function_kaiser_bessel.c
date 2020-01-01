#include "__window_functions.h"

/*  Kaiser-Bessel function with alpha = 2pi                                  */
float Kaiser_Bessel_2_0_Float(float x, float W)
{
    float bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

double Kaiser_Bessel_2_0_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

long double Kaiser_Bessel_2_0_Long_Double(long double x, long double W)
{
    long double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_0_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

/*  For all integer types, convert to double and compute.                     */
double Kaiser_Bessel_2_0_Char(char x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_UChar(unsigned char x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_Short(short x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_UShort(unsigned short x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_Int(int x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_UInt(unsigned int x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_Long(long x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_ULong(unsigned long x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_Long_Long(long long x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

double Kaiser_Bessel_2_0_ULong_Long(unsigned long long x, double W)
{
    return Kaiser_Bessel_2_0_Double((double)x, W);
}

/*  Kaiser-Bessel function with alpha = 2.5pi.                               */
float Kaiser_Bessel_2_5_Float(float x, float W)
{
    float bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

double Kaiser_Bessel_2_5_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = x*KAISER_BESSEL_2_5_A18 + KAISER_BESSEL_2_5_A17;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A16;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A15;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A14;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A13;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

long double Kaiser_Bessel_2_5_Long_Double(long double x, long double W)
{
    long double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = x*KAISER_BESSEL_2_5_A22 + KAISER_BESSEL_2_5_A21;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A20;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A19;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A18;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A17;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A16;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A15;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A14;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A13;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_2_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

/*  For all integer types, convert to double and compute.                     */
double Kaiser_Bessel_2_5_Char(char x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_UChar(unsigned char x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_Short(short x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_UShort(unsigned short x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_Int(int x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_UInt(unsigned int x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_Long(long x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_ULong(unsigned long x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_Long_Long(long long x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

double Kaiser_Bessel_2_5_ULong_Long(unsigned long long x, double W)
{
    return Kaiser_Bessel_2_5_Double((double)x, W);
}

/*  Kaiser-Bessel function with alpha = 3.5pi.                               */
float Kaiser_Bessel_3_5_Float(float x, float W)
{
    float bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = x*KAISER_BESSEL_3_5_A16 + MODIFIED_KAISER_BESSEL_3_5_A15;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A14;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A13;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

double Kaiser_Bessel_3_5_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = x*KAISER_BESSEL_3_5_A20 + MODIFIED_KAISER_BESSEL_3_5_A19;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A18;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A17;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A16;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A15;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A14;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A13;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

long double Kaiser_Bessel_3_5_Long_Double(long double x, long double W)
{
    long double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = x*KAISER_BESSEL_3_5_A24 + MODIFIED_KAISER_BESSEL_3_5_A23;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A22;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A21;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A20;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A19;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A18;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A17;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A16;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A15;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A14;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A13;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A12;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = x*bessel_x + KAISER_BESSEL_3_5_A00;
        return bessel_x;
    }
    else {
        return 0.0;
    }
}

/*  For all integer types, convert to double and compute.                     */
double Kaiser_Bessel_3_5_Char(char x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_UChar(unsigned char x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_Short(short x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_UShort(unsigned short x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_Int(int x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_UInt(unsigned int x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_Long(long x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_ULong(unsigned long x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_Long_Long(long long x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

double Kaiser_Bessel_3_5_ULong_Long(unsigned long long x, double W)
{
    return Kaiser_Bessel_3_5_Double((double)x, W);
}

/*  Kaiser-Bessel function with arbitrary alpha.                              */
float Kaiser_Bessel_Al_Float(float x, float W, float alpha)
{
    if (alpha == 0.0) return 1.0;
    else
    {
        x = 2.0*x/W;
        x = 1.0 - x*x;
        alpha *= ONE_PI;
        return BesselI0_Float(alpha*sqrtf(x))/BesselI0_Float(alpha);
    }
}

double Kaiser_Bessel_Al_Double(double x, double W, double alpha)
{
    if (alpha == 0.0) return 1.0;
    else
    {
        x = 2.0*x/W;
        x = 1.0 - x*x;
        alpha *= ONE_PI;
        return BesselI0_Double(alpha*sqrt(x))/BesselI0_Double(alpha);
    }
}

long double
Kaiser_Bessel_Al_Long_Double(long double x, long double W, long double alpha)
{
    if (alpha == 0.0) return 1.0;
    else
    {
        x = 2.0*x/W;
        x = 1.0 - x*x;
        alpha *= ONE_PI;
        return BesselI0_Long_Double(alpha*sqrtl(x))/BesselI0_Long_Double(alpha);
    }
}

/*  For all integer types, convert to double and compute.                     */
double Kaiser_Bessel_Al_Char(char x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_UChar(unsigned char x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_Short(short x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_UShort(unsigned short x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_Int(int x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_UInt(unsigned int x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_Long(long x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_ULong(unsigned long x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_Long_Long(long long x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}

double Kaiser_Bessel_Al_ULong_Long(unsigned long long x, double W, double alpha)
{
    return Kaiser_Bessel_Al_Double((double)x, W, alpha);
}
