/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/* Taylor Expansion of Kaiser-Bessel Function, alpha = 2.5 Pi                 */
#define KAISER_BESSEL_2_5_A00   2.68081686402441505171337754074e-3
#define KAISER_BESSEL_2_5_A01   4.13415654995155304033725418114e-2
#define KAISER_BESSEL_2_5_A02   1.59384725313258776340026287383e-1
#define KAISER_BESSEL_2_5_A03   2.73101421253152312766441954607e-1
#define KAISER_BESSEL_2_5_A04   2.63222948158581149121867402070e-1
#define KAISER_BESSEL_2_5_A05   1.62369147975853037148863339646e-1
#define KAISER_BESSEL_2_5_A06   6.95537872163024505198825271932e-2
#define KAISER_BESSEL_2_5_A07   2.18899350900016401453595086695e-2
#define KAISER_BESSEL_2_5_A08   5.27453612558204524028021710190e-3
#define KAISER_BESSEL_2_5_A09   1.00419724071661588837240301149e-3
#define KAISER_BESSEL_2_5_A10   1.54859836039664097565228705845e-4
#define KAISER_BESSEL_2_5_A11   1.97366389376168033180903245885e-5
#define KAISER_BESSEL_2_5_A12   2.11363735374798616835586093526e-6
#define KAISER_BESSEL_2_5_A13   1.92869494534559817370007771208e-7
#define KAISER_BESSEL_2_5_A14   1.51749490760058320452782176337e-8
#define KAISER_BESSEL_2_5_A15   1.04007461240870872184764128171e-9
#define KAISER_BESSEL_2_5_A16   6.26533506597308376631117310256e-11
#define KAISER_BESSEL_2_5_A17   3.34322980868444174824594879448e-12
#define KAISER_BESSEL_2_5_A18   1.59125943449291573576710510952e-13
#define KAISER_BESSEL_2_5_A19   6.79756800465124874234886333858e-15
#define KAISER_BESSEL_2_5_A20   2.62067605841445993538911012361e-16
#define KAISER_BESSEL_2_5_A21   9.16419924885093071618022637429e-18
#define KAISER_BESSEL_2_5_A22   2.91990641912832322368801146642e-19
#define KAISER_BESSEL_2_5_A23   8.51202777794127134835803898131e-21
#define KAISER_BESSEL_2_5_A24   2.27892650877181112915059025385e-22
#define KAISER_BESSEL_2_5_A25   5.62302577518336848756365063549e-24
#define KAISER_BESSEL_2_5_A26   1.28275332697366577347919861746e-25
#define KAISER_BESSEL_2_5_A27   2.71353478253843209761597281361e-27
#define KAISER_BESSEL_2_5_A28   5.33751491395295860622553652485e-29
#define KAISER_BESSEL_2_5_A29   9.78729947340038078323679301792e-31
#define KAISER_BESSEL_2_5_A30   1.67702732565020764195090761295e-32

float rssringoccs_Float_Kaiser_Bessel_2_5(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Kaiser_Bessel_2_5(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Kaiser_Bessel_2_5(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = KAISER_BESSEL_2_5_A12;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A11;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A10;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A09;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A08;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A07;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A06;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A05;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A04;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A03;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A02;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A01;
        bessel_x = arg*bessel_x + KAISER_BESSEL_2_5_A00;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
