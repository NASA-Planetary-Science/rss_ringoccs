/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Taylor Expansion of Modified Kaiser-Bessel Function, alpha = 3.5 Pi       */
#define MODIFIED_KAISER_BESSEL_3_5_A00  0.0
#define MODIFIED_KAISER_BESSEL_3_5_A01  4.16514995414684968723698450860e-3
#define MODIFIED_KAISER_BESSEL_3_5_A02  3.14736052127124923141946151611e-2
#define MODIFIED_KAISER_BESSEL_3_5_A03  1.05701177734383843513797151295e-1
#define MODIFIED_KAISER_BESSEL_3_5_A04  1.99680514216455319565277233641e-1
#define MODIFIED_KAISER_BESSEL_3_5_A05  2.41419041035507923786083460554e-1
#define MODIFIED_KAISER_BESSEL_3_5_A06  2.02695852544496965664760985660e-1
#define MODIFIED_KAISER_BESSEL_3_5_A07  1.25032992397232932630182925731e-1
#define MODIFIED_KAISER_BESSEL_3_5_A08  5.90500804982519735481437039452e-2
#define MODIFIED_KAISER_BESSEL_3_5_A09  2.20349118704793214846304905263e-2
#define MODIFIED_KAISER_BESSEL_3_5_A10  6.66019830971902766697459162007e-3
#define MODIFIED_KAISER_BESSEL_3_5_A11  1.66371002321116235031980438335e-3
#define MODIFIED_KAISER_BESSEL_3_5_A12  3.49213467271634973449651728719e-4
#define MODIFIED_KAISER_BESSEL_3_5_A13  6.24568564725171443222492219556e-5
#define MODIFIED_KAISER_BESSEL_3_5_A14  9.63163227374002502471798835683e-6
#define MODIFIED_KAISER_BESSEL_3_5_A15  1.29387767045843444346516264272e-6
#define MODIFIED_KAISER_BESSEL_3_5_A16  1.52766840036758485895011903343e-7
#define MODIFIED_KAISER_BESSEL_3_5_A17  1.59774363238749943364299552384e-8
#define MODIFIED_KAISER_BESSEL_3_5_A18  1.49052041225922834275380714130e-9
#define MODIFIED_KAISER_BESSEL_3_5_A19  1.24797696367142149799975692080e-10
#define MODIFIED_KAISER_BESSEL_3_5_A20  9.43023293316203461542769258304e-12
#define MODIFIED_KAISER_BESSEL_3_5_A21  6.46337975419679363032877691025e-13
#define MODIFIED_KAISER_BESSEL_3_5_A22  4.03636242526949126390645086914e-14
#define MODIFIED_KAISER_BESSEL_3_5_A23  2.30627093275607624016453712912e-15
#define MODIFIED_KAISER_BESSEL_3_5_A24  1.21021821359988898104036955063e-16
#define MODIFIED_KAISER_BESSEL_3_5_A25  5.85274375353974322653710330381e-18
#define MODIFIED_KAISER_BESSEL_3_5_A26  2.61690921774299279455177792117e-19
#define MODIFIED_KAISER_BESSEL_3_5_A27  1.08501807090336939806351306294e-20
#define MODIFIED_KAISER_BESSEL_3_5_A28  4.18308559681616420823545246848e-22
#define MODIFIED_KAISER_BESSEL_3_5_A29  1.50340710523710293537301211735e-23
#define MODIFIED_KAISER_BESSEL_3_5_A30  5.04905302598179902384690896325e-25

float rssringoccs_Float_Modified_Kaiser_Bessel_3_5(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Modified_Kaiser_Bessel_3_5(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Modified_Kaiser_Bessel_3_5(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_3_5_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_3_5_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
