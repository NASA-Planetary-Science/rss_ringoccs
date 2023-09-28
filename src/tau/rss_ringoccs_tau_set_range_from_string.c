#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_bool.h>
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
void
rssringoccs_Tau_Set_Range_From_String(const char *range,
                                      rssringoccs_TAUObj *tau)
{
    char *range_str;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (range == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Range_From_String\n\n"
            "\rInput string is NULL. Returning.\n"
        );
        return;
    }


    range_str = tmpl_String_Duplicate(range);
    tmpl_String_Remove_Whitespace(range_str);
    tmpl_String_Make_Lower_Case(range_str);

    if (strcmp(range_str, "all") == 0)
    {
        tau->rng_list[0] = 1.0;
        tau->rng_list[1] = 400000.0;
    }
    else if ((strcmp(range_str, "besselbarnard") == 0) ||
             (strcmp(range_str, "bessel-barnard") == 0))
    {
        tau->rng_list[0] = 120210.0;
        tau->rng_list[1] = 120330.0;
    }

    /*  This region is only interesting for low-inclination occultations,     *
     *  such as (Cassini) Rev133 and revs of similar occultation number.      */
    else if (strcmp(range_str, "cringripples") == 0)
    {
        tau->rng_list[0] = 77690.0;
        tau->rng_list[1] = 77760.0;
    }

    /*  Excellent test-case for most occultations since the Encke gap acts    *
     *  like a near perfect square well. Examining this region is worthwhile. */
    else if ((strcmp(range_str, "encke") == 0) ||
             (strcmp(range_str, "enckegap") == 0))
    {
        tau->rng_list[0] = 132900.0;
        tau->rng_list[1] = 134200.0;
    }
    else if ((strcmp(range_str, "herschel") == 0) ||
             (strcmp(range_str, "herschelgap") == 0))
    {
        tau->rng_list[0] = 118100.0;
        tau->rng_list[1] = 118380.0;
    }
    else if ((strcmp(range_str, "huygens") == 0) ||
             (strcmp(range_str, "huygensringlet") == 0))
    {
        tau->rng_list[0] = 117650.0;
        tau->rng_list[1] = 117950.0;
    }
    else if (strcmp(range_str, "janusepimetheus") == 0)
    {
        tau->rng_list[0] = 96200.0;
        tau->rng_list[1] = 96800.0;
    }
    else if ((strcmp(range_str, "jeffreys") == 0) ||
             (strcmp(range_str, "jeffreysgap") == 0))
    {
        tau->rng_list[0] = 118900.0;
        tau->rng_list[1] = 119000.0;
    }
    else if ((strcmp(range_str, "kuiper") == 0) ||
             (strcmp(range_str, "kuipergap") == 0))
    {
        tau->rng_list[0] = 119300.0;
        tau->rng_list[1] = 119500.0;
    }

    /*  The go-to test-case for Team-Cassini 2017-2019. An excellent gap to   *
     *  examine and test routines on. The excellent geometry of Rev007 means  *
     *  one can get very accurate results on the Maxwell ringlet using any of *
     *  the reconstruction methods.                                           */
    else if ((strcmp(range_str, "maxwell") == 0) ||
             (strcmp(range_str, "maxwellringlet") == 0))
    {
        tau->rng_list[0] = 87410.0;
        tau->rng_list[1] = 87610.0;
    }
    else if ((strcmp(range_str, "russell") == 0) ||
             (strcmp(range_str, "russellgap") == 0))
    {
        tau->rng_list[0] = 118550.0;
        tau->rng_list[1] = 118660.0;
    }

    /*  A very steep gap, excellent for testing the resolution of the         *
     *  reconstruction.                                                       */
    else if ((strcmp(range_str, "titan") == 0) ||
             (strcmp(range_str, "titanringlet") == 0))
    {
        tau->rng_list[0] = 77870.0;
        tau->rng_list[1] = 77930.0;
    }
    else
    {
        char errmes1[1024];
        const char *errmes2 =
            "\r\thuygensringlet    [117650.0, 117950.0]\n"
            "\r\tjanusepimetheus   [96200.0, 96800.0]\n"
            "\r\tjeffreys          [118900.0, 119000.0]\n"
            "\r\tjeffreysgap       [118900.0, 119000.0]\n"
            "\r\tkuiper            [119300.0, 119500.0]\n"
            "\r\tkuipergap         [119300.0, 119500.0]\n"
            "\r\tmaxwell           [87410.0, 87610.0]\n"
            "\r\tmaxwellringlet    [87410.0, 87610.0]\n"
            "\r\trussell           [118550.0, 118660.0]\n"
            "\r\trussellgap        [118550.0, 118660.0]\n"
            "\r\ttitan             [77870.0, 77930.0]\n"
            "\r\ttitanringlet      [77870.0, 77930.0\n\n";

        strcpy(
            errmes1,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Range_From_String\n\n"
            "\rIllegal string for range. Allowed strings are:\n"
            "\r\tall               [1.0, 400000.0]\n"
            "\r\tbesselbarnard     [120210.0, 120330.0]\n"
            "\r\tbessel-barnard    [120210.0, 120330.0]\n"
            "\r\tcringripples      [77690.0, 77760.0]\n"
            "\r\tencke             [132900.0, 134200.0]\n"
            "\r\tenckegap          [132900.0, 134200.0]\n"
            "\r\therschel          [118100.0, 118380.0]\n"
            "\r\therschelgap       [118100.0, 118380.0]\n"
            "\r\thuygens           [117650.0, 117950.0]\n"
        );

        strcat(errmes1, errmes2);

        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(errmes1);
        tau->rng_list[0] = -1.0;
        tau->rng_list[1] = -1.0;
    }
}
