/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Taylor Expansion of Modified Kaiser-Bessel Function, alpha = 2.5 Pi       */
#define MODIFIED_KAISER_BESSEL_2_5_A00  0.0
#define MODIFIED_KAISER_BESSEL_2_5_A01  4.14526925768347286292556360856e-2
#define MODIFIED_KAISER_BESSEL_2_5_A02  1.59813155114582899876788606369e-1
#define MODIFIED_KAISER_BESSEL_2_5_A03  2.73835524144247175471513296242e-1
#define MODIFIED_KAISER_BESSEL_2_5_A04  2.63930497487175121688556222589e-1
#define MODIFIED_KAISER_BESSEL_2_5_A05  1.62805599973820467924927444509e-1
#define MODIFIED_KAISER_BESSEL_2_5_A06  6.97407493933859400420106113296e-2
#define MODIFIED_KAISER_BESSEL_2_5_A07  2.19487757381451495909574422227e-2
#define MODIFIED_KAISER_BESSEL_2_5_A08  5.28871419979786875646449115789e-3
#define MODIFIED_KAISER_BESSEL_2_5_A09  1.00689654595734627672697987971e-3
#define MODIFIED_KAISER_BESSEL_2_5_A10  1.55276102834723413370648946123e-4
#define MODIFIED_KAISER_BESSEL_2_5_A11  1.97896914762602018893080108449e-5
#define MODIFIED_KAISER_BESSEL_2_5_A12  2.11931885948674321924199402729e-6
#define MODIFIED_KAISER_BESSEL_2_5_A13  1.93387932164404963793836133629e-7
#define MODIFIED_KAISER_BESSEL_2_5_A14  1.52157396875588453820420936413e-8
#define MODIFIED_KAISER_BESSEL_2_5_A15  1.04287035684833890841740681814e-9
#define MODIFIED_KAISER_BESSEL_2_5_A16  6.28217643049071999192344641396e-11
#define MODIFIED_KAISER_BESSEL_2_5_A17  3.35221648717512149821851074143e-12
#define MODIFIED_KAISER_BESSEL_2_5_A18  1.59553677639023382676771341351e-13
#define MODIFIED_KAISER_BESSEL_2_5_A19  6.81584002352881736180902151748e-15
#define MODIFIED_KAISER_BESSEL_2_5_A20  2.62772049583363333730995733786e-16
#define MODIFIED_KAISER_BESSEL_2_5_A21  9.18883282685386168648152273945e-18
#define MODIFIED_KAISER_BESSEL_2_5_A22  2.92775519462781655781417909659e-19
#define MODIFIED_KAISER_BESSEL_2_5_A23  8.53490830405568585080958796776e-21
#define MODIFIED_KAISER_BESSEL_2_5_A24  2.28505231555452774631491924009e-22
#define MODIFIED_KAISER_BESSEL_2_5_A25  5.63814059757909922984721163965e-24
#define MODIFIED_KAISER_BESSEL_2_5_A26  1.28620139737026775520317332179e-25
#define MODIFIED_KAISER_BESSEL_2_5_A27  2.72082882634020868563820934433e-27
#define MODIFIED_KAISER_BESSEL_2_5_A28  5.35186227659799861946749973142e-29
#define MODIFIED_KAISER_BESSEL_2_5_A29  9.81360795911409807961459041065e-31
#define MODIFIED_KAISER_BESSEL_2_5_A30  1.68153521360829966876248982166e-32

float rssringoccs_Float_Modified_Kaiser_Bessel_2_5(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Modified_Kaiser_Bessel_2_5(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Modified_Kaiser_Bessel_2_5(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
