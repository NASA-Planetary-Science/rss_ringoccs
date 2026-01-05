#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Error message for the legal psitypes. Defined at the bottom of this file. */
static const char rssringoccs_range_string_error_message[969];

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
        tau->error_message =
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Range_From_String\n\n"
            "\rInput string is NULL.\n\n";

        return;
    }


    range_str = tmpl_String_Duplicate(range);
    tmpl_String_Remove_Whitespace(range_str);
    tmpl_String_Make_Lower_Case(range_str);

    if (tmpl_String_Are_Equal(range_str, "all"))
    {
        tau->rng_list[0] = 1.0;
        tau->rng_list[1] = 400000.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "besselbarnard")) ||
             (tmpl_String_Are_Equal(range_str, "bessel-barnard")))
    {
        tau->rng_list[0] = 120210.0;
        tau->rng_list[1] = 120330.0;
    }

    /*  This region is only interesting for low-inclination occultations,     *
     *  such as (Cassini) Rev133 and revs of similar occultation number.      */
    else if (tmpl_String_Are_Equal(range_str, "cringripples"))
    {
        tau->rng_list[0] = 77690.0;
        tau->rng_list[1] = 77760.0;
    }

    /*  Excellent test-case for most occultations since the Encke gap acts    *
     *  like a near perfect square well. Examining this region is worthwhile. */
    else if ((tmpl_String_Are_Equal(range_str, "encke")) ||
             (tmpl_String_Are_Equal(range_str, "enckegap")))
    {
        tau->rng_list[0] = 132900.0;
        tau->rng_list[1] = 134200.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "herschel")) ||
             (tmpl_String_Are_Equal(range_str, "herschelgap")))
    {
        tau->rng_list[0] = 118100.0;
        tau->rng_list[1] = 118380.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "huygens")) ||
             (tmpl_String_Are_Equal(range_str, "huygensringlet")))
    {
        tau->rng_list[0] = 117650.0;
        tau->rng_list[1] = 117950.0;
    }
    else if (tmpl_String_Are_Equal(range_str, "janusepimetheus"))
    {
        tau->rng_list[0] = 96200.0;
        tau->rng_list[1] = 96800.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "jeffreys")) ||
             (tmpl_String_Are_Equal(range_str, "jeffreysgap")))
    {
        tau->rng_list[0] = 118900.0;
        tau->rng_list[1] = 119000.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "kuiper")) ||
             (tmpl_String_Are_Equal(range_str, "kuipergap")))
    {
        tau->rng_list[0] = 119300.0;
        tau->rng_list[1] = 119500.0;
    }

    /*  The go-to test-case for Team-Cassini 2017-2019. An excellent gap to   *
     *  examine and test routines on. The excellent geometry of Rev007 means  *
     *  one can get very accurate results on the Maxwell ringlet using any of *
     *  the reconstruction methods.                                           */
    else if ((tmpl_String_Are_Equal(range_str, "maxwell")) ||
             (tmpl_String_Are_Equal(range_str, "maxwellringlet")))
    {
        tau->rng_list[0] = 87410.0;
        tau->rng_list[1] = 87610.0;
    }
    else if ((tmpl_String_Are_Equal(range_str, "russell")) ||
             (tmpl_String_Are_Equal(range_str, "russellgap")))
    {
        tau->rng_list[0] = 118550.0;
        tau->rng_list[1] = 118660.0;
    }

    /*  A very steep gap, excellent for testing the resolution of the         *
     *  reconstruction.                                                       */
    else if ((tmpl_String_Are_Equal(range_str, "titan")) ||
             (tmpl_String_Are_Equal(range_str, "titanringlet")))
    {
        tau->rng_list[0] = 77870.0;
        tau->rng_list[1] = 77930.0;
    }
    else
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = rssringoccs_range_string_error_message;
        tau->rng_list[0] = -1.0;
        tau->rng_list[1] = -1.0;
    }
}

/*  The string for the error message. C89 compilers are required to support   *
 *  string literals of up to ~500 characters. The range error message is      *
 *  about 1000 characters long, so we provide it as a char array. This makes  *
 *  it infinitely less readable, but portable. The array below was created    *
 *  using a short C program, it was not typed out by hand.                    */
static const char rssringoccs_range_string_error_message[969] = {
    '\n', '\r', 'E', 'r', 'r', 'o', 'r', ' ', 'E', 'n', 'c', 'o', 'u', 'n',
    't', 'e', 'r', 'e', 'd', ':', ' ', 'r', 's', 's', '_', 'r', 'i', 'n', 'g',
    'o', 'c', 'c', 's', '\n', '\r', '\t', 'r', 's', 's', 'r', 'i', 'n', 'g',
    'o', 'c', 'c', 's', '_', 'T', 'a', 'u', '_', 'G', 'e', 't', '_', 'R', 'a',
    'n', 'g', 'e', '_', 'F', 'r', 'o', 'm', '_', 'S', 't', 'r', 'i', 'n', 'g',
    '\n', '\n', '\r', 'I', 'l', 'l', 'e', 'g', 'a', 'l', ' ', 's', 't', 'r',
    'i', 'n', 'g', ' ', 'f', 'o', 'r', ' ', 'r', 'a', 'n', 'g', 'e', '.', ' ',
    'A', 'l', 'l', 'o', 'w', 'e', 'd', ' ', 's', 't', 'r', 'i', 'n', 'g', 's',
    ' ', 'a', 'r', 'e', ':', '\n', '\r', '\t', 'a', 'l', 'l', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '.',
    '0', ',', ' ', '4', '0', '0', '0', '0', '0', '.', '0', ']', '\n', '\r',
    '\t', 'b', 'e', 's', 's', 'e', 'l', 'b', 'a', 'r', 'n', 'a', 'r', 'd',
    ' ', ' ', ' ', ' ', ' ', '[', '1', '2', '0', '2', '1', '0', '.', '0', ',',
    ' ', '1', '2', '0', '3', '3', '0', '.', '0', ']', '\n', '\r', '\t', 'b',
    'e', 's', 's', 'e', 'l', '-', 'b', 'a', 'r', 'n', 'a', 'r', 'd', ' ', ' ',
    ' ', ' ', '[', '1', '2', '0', '2', '1', '0', '.', '0', ',', ' ', '1', '2',
    '0', '3', '3', '0', '.', '0', ']', '\n', '\r', '\t', 'c', 'r', 'i', 'n',
    'g', 'r', 'i', 'p', 'p', 'l', 'e', 's', ' ', ' ', ' ', ' ', ' ', ' ', '[',
    '7', '7', '6', '9', '0', '.', '0', ',', ' ', '7', '7', '7', '6', '0', '.',
    '0', ']', '\n', '\r', '\t', 'e', 'n', 'c', 'k', 'e', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '3', '2', '9', '0',
    '0', '.', '0', ',', ' ', '1', '3', '4', '2', '0', '0', '.', '0', ']', '\n',
    '\r', '\t', 'e', 'n', 'c', 'k', 'e', 'g', 'a', 'p', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '3', '2', '9', '0', '0', '.', '0',
    ',', ' ', '1', '3', '4', '2', '0', '0', '.', '0', ']', '\n', '\r', '\t',
    'h', 'e', 'r', 's', 'c', 'h', 'e', 'l', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', '[', '1', '1', '8', '1', '0', '0', '.', '0', ',', ' ', '1',
    '1', '8', '3', '8', '0', '.', '0', ']', '\n', '\r', '\t', 'h', 'e', 'r',
    's', 'c', 'h', 'e', 'l', 'g', 'a', 'p', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    '[', '1', '1', '8', '1', '0', '0', '.', '0', ',', ' ', '1', '1', '8', '3',
    '8', '0', '.', '0', ']', '\n', '\r', '\t', 'h', 'u', 'y', 'g', 'e', 'n',
    's', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '1',
    '7', '6', '5', '0', '.', '0', ',', ' ', '1', '1', '7', '9', '5', '0', '.',
    '0', ']', '\n', '\r', '\t', 'h', 'u', 'y', 'g', 'e', 'n', 's', 'r', 'i',
    'n', 'g', 'l', 'e', 't', ' ', ' ', ' ', ' ', '[', '1', '1', '7', '6', '5',
    '0', '.', '0', ',', ' ', '1', '1', '7', '9', '5', '0', '.', '0', ']', '\n',
    '\r', '\t', 'j', 'a', 'n', 'u', 's', 'e', 'p', 'i', 'm', 'e', 't', 'h',
    'e', 'u', 's', ' ', ' ', ' ', '[', '9', '6', '2', '0', '0', '.', '0', ',',
    ' ', '9', '6', '8', '0', '0', '.', '0', ']', '\n', '\r', '\t', 'j', 'e',
    'f', 'f', 'r', 'e', 'y', 's', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', '[', '1', '1', '8', '9', '0', '0', '.', '0', ',', ' ', '1', '1', '9',
    '0', '0', '0', '.', '0', ']', '\n', '\r', '\t', 'j', 'e', 'f', 'f', 'r',
    'e', 'y', 's', 'g', 'a', 'p', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1',
    '1', '8', '9', '0', '0', '.', '0', ',', ' ', '1', '1', '9', '0', '0', '0',
    '.', '0', ']', '\n', '\r', '\t', 'k', 'u', 'i', 'p', 'e', 'r', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '1', '9', '3',
    '0', '0', '.', '0', ',', ' ', '1', '1', '9', '5', '0', '0', '.', '0', ']',
    '\n', '\r', '\t', 'k', 'u', 'i', 'p', 'e', 'r', 'g', 'a', 'p', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '1', '9', '3', '0', '0', '.',
    '0', ',', ' ', '1', '1', '9', '5', '0', '0', '.', '0', ']', '\n', '\r',
    '\t', 'm', 'a', 'x', 'w', 'e', 'l', 'l', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', '[', '8', '7', '4', '1', '0', '.', '0', ',', ' ',
    '8', '7', '6', '1', '0', '.', '0', ']', '\n', '\r', '\t', 'm', 'a', 'x',
    'w', 'e', 'l', 'l', 'r', 'i', 'n', 'g', 'l', 'e', 't', ' ', ' ', ' ', ' ',
    '[', '8', '7', '4', '1', '0', '.', '0', ',', ' ', '8', '7', '6', '1', '0',
    '.', '0', ']', '\n', '\r', '\t', 'r', 'u', 's', 's', 'e', 'l', 'l', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '1', '8', '5',
    '5', '0', '.', '0', ',', ' ', '1', '1', '8', '6', '6', '0', '.', '0', ']',
    '\n', '\r', '\t', 'r', 'u', 's', 's', 'e', 'l', 'l', 'g', 'a', 'p', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', '[', '1', '1', '8', '5', '5', '0', '.',
    '0', ',', ' ', '1', '1', '8', '6', '6', '0', '.', '0', ']', '\n', '\r',
    '\t', 't', 'i', 't', 'a', 'n', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', '[', '7', '7', '8', '7', '0', '.', '0', ',', ' ',
    '7', '7', '9', '3', '0', '.', '0', ']', '\n', '\r', '\t', 't', 'i', 't',
    'a', 'n', 'r', 'i', 'n', 'g', 'l', 'e', 't', ' ', ' ', ' ', ' ', ' ', ' ',
    '[', '7', '7', '8', '7', '0', '.', '0', ',', ' ', '7', '7', '9', '3', '0',
    '.', '0', ']', '\n', '\n', '\0'
};
