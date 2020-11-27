/*  The speed of light is defined here.                                       */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Convert a frequency (in hertz) to a wavelength (in kilometers).           */
#define _define_freq_to_wave(type, Type)                                       \
type rssringoccs_##Type##_Frequency_To_Wavelength(type frequency)              \
{                                                                              \
    /*  The conversion is the speed of light divided by the input frequency. */\
    return SPEED_OF_LIGHT_KMS / frequency;                                     \
}

_define_freq_to_wave(float, Float)
_define_freq_to_wave(double, Double)
_define_freq_to_wave(long double, LongDouble)
