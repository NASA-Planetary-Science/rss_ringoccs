#ifndef RSS_RINGOCCS_MODIFIED_KAISER_BESSEL_H
#define RSS_RINGOCCS_MODIFIED_KAISER_BESSEL_H

#include "__math_constants.h"

/*  Modified Kaiser-Bessel with alpha = 2.0pi                                */
float Modified_Kaiser_Bessel_2_0_Float(float x, float W)
{
    float bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

double Modified_Kaiser_Bessel_2_0_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

long double Modified_Kaiser_Bessel_2_0_Long_Double(long double x, long double W)
{
    long double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

/*  Modified Kaiser-Bessel with alpha = 2.5pi                                */
float Modified_Kaiser_Bessel_2_5_Float(float x, float W)
{
    float bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

double Modified_Kaiser_Bessel_2_5_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

long double Modified_Kaiser_Bessel_2_5_Long_Double(long double x,long double W)
{
    long double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x > 0){
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

/*  Modified Kaiser-Bessel with alpha = 2.5pi                                */
double Modified_Kaiser_Bessel_3_5_Double(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = MODIFIED_KAISER_BESSEL_3_5_A22;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A21;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A20;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A19;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A18;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A17;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A16;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A15;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A14;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A13;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A12;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A11;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A10;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A09;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A08;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A07;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A06;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A05;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A04;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A03;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A02;
        bessel_x = x*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A01;
        return bessel_x*x;
    }
    else {
        return 0.0;
    }
}

#endif