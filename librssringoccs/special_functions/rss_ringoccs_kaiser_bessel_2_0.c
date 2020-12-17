/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/* Taylor Expansion of Kaiser-Bessel Function, alpha = 2.0 Pi                 */
#define KAISER_BESSEL_2_0_A00   1.14799345379586479212325625173e-2
#define KAISER_BESSEL_2_0_A01   1.13302412440054404174028755736e-1
#define KAISER_BESSEL_2_0_A02   2.79562497118100660800421522589e-1
#define KAISER_BESSEL_2_0_A03   3.06574583548481935246360365287e-1
#define KAISER_BESSEL_2_0_A04   1.89110616190764661645757930810e-1
#define KAISER_BESSEL_2_0_A05   7.46578787939636569441597855450e-2
#define KAISER_BESSEL_2_0_A06   2.04678813644694334928965293809e-2
#define KAISER_BESSEL_2_0_A07   4.12265085705596709062582845471e-3
#define KAISER_BESSEL_2_0_A08   6.35764578796162329491799908822e-4
#define KAISER_BESSEL_2_0_A09   7.74659862338682228855283637288e-5
#define KAISER_BESSEL_2_0_A10   7.64558638668513481381569815368e-6
#define KAISER_BESSEL_2_0_A11   6.23627380586252011565352895678e-7
#define KAISER_BESSEL_2_0_A12   4.27427468060687589919240887794e-8
#define KAISER_BESSEL_2_0_A13   2.49617752657884230813659743142e-9
#define KAISER_BESSEL_2_0_A14   1.25695330113382099459238476608e-10
#define KAISER_BESSEL_2_0_A15   5.51361414792629226182398674779e-12
#define KAISER_BESSEL_2_0_A16   2.12567150235476125296700259059e-13
#define KAISER_BESSEL_2_0_A17   7.25935529927708671111891419534e-15
#define KAISER_BESSEL_2_0_A18   2.21132608058075595943487946625e-16
#define KAISER_BESSEL_2_0_A19   6.04568244242202512396089987671e-18
#define KAISER_BESSEL_2_0_A20   1.49171235103292705239359580748e-19
#define KAISER_BESSEL_2_0_A21   3.33846049544533584752889898662e-21
#define KAISER_BESSEL_2_0_A22   6.80770338816327424681834548673e-23
#define KAISER_BESSEL_2_0_A23   1.27011983593813213204479126306e-24
#define KAISER_BESSEL_2_0_A24   2.17631602824407691118793520341e-26
#define KAISER_BESSEL_2_0_A25   3.43670052008304870068386995834e-28
#define KAISER_BESSEL_2_0_A26   5.01758499678073153159604331079e-30
#define KAISER_BESSEL_2_0_A27   6.79308353457709783558803425981e-32
#define KAISER_BESSEL_2_0_A28   8.55166417727420676456688381397e-34
#define KAISER_BESSEL_2_0_A29   1.00358552200551391399489779147e-35
#define KAISER_BESSEL_2_0_A30   1.10055467609502016532883145908e-37

float rssringoccs_Float_Kaiser_Bessel_2_0(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Kaiser_Bessel_2_0(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Kaiser_Bessel_2_0(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_0_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
