#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

/*  Various coefficients and constants defined here.    */
#include "__math_constants.h"

double __rect(double x, double W)
{
    if (fabs(x) <= W/2.0){
        return 1.0;
    }
    else {
        return 0.0;
    }
}

double __coss(double x, double W)
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

double Kaiser_Bessel_Window_2_0(double x, double W)
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

double Kaiser_Bessel_Window_2_5(double x, double W)
{
    double bessel_x;
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

double Kaiser_Bessel_Window_3_5(double x, double W)
{
    double bessel_x;
    x = 2.0*x/W;
    x = 1.0 - x*x;
    if (x >= 0){
        bessel_x = KAISER_BESSEL_3_5_A22;
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

double Modified_Kaiser_Bessel_Window_2_0(double x, double W)
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

double Modified_Kaiser_Bessel_Window_2_5(double x, double W)
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

double Modified_Kaiser_Bessel_Window_3_5(double x, double W)
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

/*
    static void __coss(double* wfunc, double w_width, double dx, long nw_pts)
    {
            long i;
            double x;
            dx = ONE_PI * dx / w_width;

            for (i=0; i<nw_pts; i++){
                x = (i-nw_pts) * dx;
                x = cos(x);
                x *= x;
                wfunc[i] = x;
            }
    }

    static void __kbmd20(double* wfunc, double w_width, double dx, long nw_pts)
    {
        long i;
        double x;

        for (i=0; i<nw_pts; i++){
            x = (i-nw_pts) * dx;
            wfunc[i] = Modified_Kaiser_Bessel_Window_2_0(x, w_width);
        }
    }
*/

#endif