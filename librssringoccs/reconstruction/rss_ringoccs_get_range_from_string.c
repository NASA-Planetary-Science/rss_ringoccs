#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  This function provides a plethora of strings for the user to provide for  *
 *  the sake of convenience in specifying features of the rings of Saturn.    *
 *  This takes in a string and a pointer to the range variable, and on        *
 *  success it changes the values rng_list points to and then returns 1. On   *
 *  failure it return 0. Hence you may use this function like:                *
 *      if (!GetRangeFromString(my_string, my_range))                         *
 *          RaiseAnError;                                                     *
 *  All of the features correspond to Saturn. The user can append their own   *
 *  features for other planets, if desired, but MAKE SURE YOU ALTER THE ERROR *
 *  CHECKING ROUTINE found in the Diffrec_init function. If you do not, even  *
 *  if you've added a new string to this function, Diffrec_init will refuse   *
 *  to accept your new string and raise an error.                             */
void GetRangeFromString(char *range, double *rng_list)
{
    if (strcmp(range, "all") == 0)
    {
        rng_list[0] = 1.0;
        rng_list[1] = 400000.0;
    }
    else if ((strcmp(range, "besselbarnard") == 0) ||
             (strcmp(range, "bessel-barnard") == 0))
    {
        rng_list[0] = 120210.0;
        rng_list[1] = 120330.0;
    }

    /*  This region is only interesting for low-inclination occultations,     *
     *  such as (Cassini) Rev133 and revs of similar occultation number.      */
    else if (strcmp(range, "cringripples") == 0)
    {
        rng_list[0] = 77690.0;
        rng_list[1] = 77760.0;
    }

    /*  Excellent test-case for most occultations since the Encke gap acts    *
     *  like a near perfect square well. Examining this region is worthwhile. */
    else if ((strcmp(range, "encke") == 0) || (strcmp(range, "enckegap") == 0))
    {
        rng_list[0] = 132900.0;
        rng_list[1] = 134200.0;
    }
    else if ((strcmp(range, "herschel") == 0) ||
             (strcmp(range, "herschelgap") == 0))
    {
        rng_list[0] = 118100.0;
        rng_list[1] = 118380.0;
    }
    else if ((strcmp(range, "huygens") == 0) ||
             (strcmp(range, "huygensringlet") == 0))
    {
        rng_list[0] = 117650.0;
        rng_list[1] = 117950.0;
    }
    else if (strcmp(range, "janusepimetheus") == 0)
    {
        rng_list[0] = 96200.0;
        rng_list[1] = 96800.0;
    }
    else if ((strcmp(range, "jeffreys") == 0) ||
             (strcmp(range, "jeffreysgap") == 0))
    {
        rng_list[0] = 118900.0;
        rng_list[1] = 119000.0;
    }
    else if ((strcmp(range, "kuiper") == 0) ||
             (strcmp(range, "kuipergap") == 0))
    {
        rng_list[0] = 119300.0;
        rng_list[1] = 119500.0;
    }

    /*  The go-to test-case for Team-Cassini 2017-2019. An excellent gap to   *
     *  examine and test routines on. The excellent geometry of Rev007 means  *
     *  one can get very accurate results on the Maxwell ringlet using any of *
     *  the reconstruction methods.                                           */
    else if ((strcmp(range, "maxwell") == 0) ||
             (strcmp(range, "maxwellringlet") == 0))
    {
        rng_list[0] = 87410.0;
        rng_list[1] = 87610.0;
    }
    else if ((strcmp(range, "russell") == 0) ||
             (strcmp(range, "russellgap") == 0))
    {
        rng_list[0] = 118550.0;
        rng_list[1] = 118660.0;
    }

    /*  A very steep gap, excellent for testing the resolution of the         *
     *  reconstruction.                                                       */
    else if ((strcmp(range, "titan") == 0) ||
             (strcmp(range, "titanringlet") == 0))
    {
        rng_list[0] = 77870.0;
        rng_list[1] = 77930.0;
    }
    else
    {
        rng_list[0] = -1.0;
        rng_list[1] = -1.0;
    }
}