/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/* Taylor Expansion of Modified Kaiser-Bessel Function, alpha = 2.0 Pi        */
#define MODIFIED_KAISER_BESSEL_2_0_A00 0.0
#define MODIFIED_KAISER_BESSEL_2_0_A01 1.14618222126928749645196222695e-1
#define MODIFIED_KAISER_BESSEL_2_0_A02 2.82809127387243423533545769115e-1
#define MODIFIED_KAISER_BESSEL_2_0_A03 3.10134912036597639686900220959e-1
#define MODIFIED_KAISER_BESSEL_2_0_A04 1.91306805797991572819499337886e-1
#define MODIFIED_KAISER_BESSEL_2_0_A05 7.55248996984881942399995746237e-2
#define MODIFIED_KAISER_BESSEL_2_0_A06 2.07055800682230959900642506315e-2
#define MODIFIED_KAISER_BESSEL_2_0_A07 4.17052824833556674429483596109e-3
#define MODIFIED_KAISER_BESSEL_2_0_A08 6.43147874291253173568797655869e-4
#define MODIFIED_KAISER_BESSEL_2_0_A09 7.83656184031631837458373758816e-5
#define MODIFIED_KAISER_BESSEL_2_0_A10 7.73437652285948594043705662616e-6
#define MODIFIED_KAISER_BESSEL_2_0_A11 6.30869723716497460710444914139e-7
#define MODIFIED_KAISER_BESSEL_2_0_A12 4.32391291785164657063120476020e-8
#define MODIFIED_KAISER_BESSEL_2_0_A13 2.52516627005655274230166467041e-9
#define MODIFIED_KAISER_BESSEL_2_0_A14 1.27155061900166085890453190618e-10
#define MODIFIED_KAISER_BESSEL_2_0_A15 5.57764514911408455397588384071e-12
#define MODIFIED_KAISER_BESSEL_2_0_A16 2.15035746528949534499532181249e-13
#define MODIFIED_KAISER_BESSEL_2_0_A17 7.34366003575659549104967048795e-15
#define MODIFIED_KAISER_BESSEL_2_0_A18 2.23700677188294233151931611484e-16
#define MODIFIED_KAISER_BESSEL_2_0_A19 6.11589248782343117103908763044e-18
#define MODIFIED_KAISER_BESSEL_2_0_A20 1.50903598536028707255650259738e-19
#define MODIFIED_KAISER_BESSEL_2_0_A21 3.37723088492383361300462693638e-21
#define MODIFIED_KAISER_BESSEL_2_0_A22 6.88676297630974740993061496477e-23
#define MODIFIED_KAISER_BESSEL_2_0_A23 1.28487006011807057474580521836e-24
#define MODIFIED_KAISER_BESSEL_2_0_A24 2.20159013891827398230794077153e-26
#define MODIFIED_KAISER_BESSEL_2_0_A25 3.47661179591403669266617928165e-28
#define MODIFIED_KAISER_BESSEL_2_0_A26 5.07585548547815911063908063335e-30
#define MODIFIED_KAISER_BESSEL_2_0_A27 6.87197333865141319387208768561e-32
#define MODIFIED_KAISER_BESSEL_2_0_A28 8.65097682491361238639269885487e-34
#define MODIFIED_KAISER_BESSEL_2_0_A29 1.01524041551592668958612323193e-35
#define MODIFIED_KAISER_BESSEL_2_0_A30 1.11333569701553102471644674153e-37

float rssringoccs_Float_Modified_Kaiser_Bessel_2_0(float x, float W)
{
    float bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

double rssringoccs_Double_Modified_Kaiser_Bessel_2_0(double x, double W)
{
    double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    /*  arg > 0 means x is within the window. We use a Taylor series for the  *
     *  modified Kaiser-Bessel function.                                      */
    if (arg > 0.0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}

long double
rssringoccs_LDouble_Modified_Kaiser_Bessel_2_0(long double x, long double W)
{
    long double bessel_x, arg;
    arg = 2.0*x/W;
    arg = 1.0 - arg*arg;

    if (arg > 0)
    {
        bessel_x = MODIFIED_KAISER_BESSEL_2_0_A12;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A11;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A10;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A09;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A08;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A07;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A06;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A05;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A04;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A03;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A02;
        bessel_x = arg*bessel_x + MODIFIED_KAISER_BESSEL_2_0_A01;
        bessel_x = arg*bessel_x;
    }

    /*  For values outside of the window, we return 0.                        */
    else
        bessel_x = 0.0;

    return bessel_x;
}
