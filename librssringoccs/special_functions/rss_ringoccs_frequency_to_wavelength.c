/*  The speed of light is defined here.                                       */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Convert a frequency (in hertz) to a wavelength (in kilometers).           */
float rssringoccs_Float_Frequency_To_Wavelength(float frequency)
{
    /*  The conversion is the speed of light divided by the input frequency.  */
    return rssringoccs_Speed_Of_Light_KMS_F / frequency;
}

double rssringoccs_Double_Frequency_To_Wavelength(double frequency)
{
    /*  The conversion is the speed of light divided by the input frequency.  */
    return rssringoccs_Speed_Of_Light_KMS / frequency;
}

long double rssringoccs_LDouble_Frequency_To_Wavelength(long double frequency)
{
    /*  The conversion is the speed of light divided by the input frequency.  */
    return rssringoccs_Speed_Of_Light_KMS_L / frequency;
}
