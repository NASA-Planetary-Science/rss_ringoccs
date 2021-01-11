/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/* Taylor Expansion of Kaiser-Bessel Function, alpha = 3.5 Pi                 */
#define KAISER_BESSEL_3_5_A00   1.37782784195108944072394913183e-4
#define KAISER_BESSEL_3_5_A01   4.16457606818957720393279069168e-3
#define KAISER_BESSEL_3_5_A02   3.14692686917576270929327830656e-2
#define KAISER_BESSEL_3_5_A03   1.05686613931822898050138235809e-1
#define KAISER_BESSEL_3_5_A04   1.99653001679257065318455314636e-1
#define KAISER_BESSEL_3_5_A05   2.41385777647876338247314746026e-1
#define KAISER_BESSEL_3_5_A06   2.02667924545588583615194331816e-1
#define KAISER_BESSEL_3_5_A07   1.25015765003424195987707818786e-1
#define KAISER_BESSEL_3_5_A08   5.90419444137539790852413424858e-2
#define KAISER_BESSEL_3_5_A09   2.20318758389723129879141296234e-2
#define KAISER_BESSEL_3_5_A10   6.65928064905262302082534080707e-3
#define KAISER_BESSEL_3_5_A11   1.66348079261207100704542994789e-3
#define KAISER_BESSEL_3_5_A12   3.49165351667835860027922841080e-4
#define KAISER_BESSEL_3_5_A13   6.24482509929402865993003387209e-5
#define KAISER_BESSEL_3_5_A14   9.63030520062900565637523002718e-6
#define KAISER_BESSEL_3_5_A15   1.29369939639059079870521852670e-6
#define KAISER_BESSEL_3_5_A16   1.52745791396205532471640815459e-7
#define KAISER_BESSEL_3_5_A17   1.59752349082139907733121625093e-8
#define KAISER_BESSEL_3_5_A18   1.49031504420692762467929656674e-9
#define KAISER_BESSEL_3_5_A19   1.24780501393075549121691645926e-10
#define KAISER_BESSEL_3_5_A20   9.42893360941289514158532861815e-12
#define KAISER_BESSEL_3_5_A21   6.46248921173895009721953651274e-13
#define KAISER_BESSEL_3_5_A22   4.03580628401651711103915833816e-14
#define KAISER_BESSEL_3_5_A23   2.30595316832585285710081611124e-15
#define KAISER_BESSEL_3_5_A24   1.21005146636493555727671823087e-16
#define KAISER_BESSEL_3_5_A25   5.85193734621019998820702643460e-18
#define KAISER_BESSEL_3_5_A26   2.61654865270498632041384548487e-19
#define KAISER_BESSEL_3_5_A27   1.08486857409265832570275404458e-20
#define KAISER_BESSEL_3_5_A28   4.18250923963610841818088138643e-22
#define KAISER_BESSEL_3_5_A29   1.50319996162036465817766810461e-23
#define KAISER_BESSEL_3_5_A30   5.04835735339833051183961001462e-25

float rssringoccs_Float_Kaiser_Bessel_3_5(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Kaiser_Bessel_3_5(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Kaiser_Bessel_3_5(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_3_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
