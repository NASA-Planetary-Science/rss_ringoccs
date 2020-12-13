/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Convert a wavelength to a wavelength.                                     */
#define _define_wavelength_to_wavenumber(type, Type)                           \
type rssringoccs_##Type##_Wavelength_To_Wavenumber(type lambda)                \
{                                                                              \
    /*  The conversion is the two pi over the wavelength.                    */\
    return rssringoccs_Two_Pi / lambda;                                        \
}

_define_wavelength_to_wavenumber(float, Float)
_define_wavelength_to_wavenumber(double, Double)
_define_wavelength_to_wavenumber(long double, LDouble)
