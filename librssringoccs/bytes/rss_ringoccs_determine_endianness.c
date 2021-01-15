/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *                   rss_ringoccs_determine_endianness                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains code for determining the endianness of a platform.           *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Determine_Endianness:                                     *
 *  Purpose:                                                                  *
 *      Determines the endianness of a platform. That is, whether the system  *
 *      is little-endian, big-endian, or mixed-endian.                        *
 *  Arguments:                                                                *
 *      None (void).                                                          *
 *  Output:                                                                   *
 *      endianness (rssringoccs_Endian):                                      *
 *          An enum data-type whose value corresponds to the endianness.      *
 *  Method:                                                                   *
 *      Use a union data type between an unsigned int and a char array whose  *
 *      size makes it so that the entire array has the same number of bits as *
 *      an unsigned int (determined via macros in limits.h). Set the          *
 *      unsigned bit portion to 1 and then see which element of the char      *
 *      array is one. This will tell us the endianness.                       *
 *  NOTES:                                                                    *
 *      You could define the function without the macros from limits.h and    *
 *      instead use sizeof. This will give compiler warnings with gcc and the *
 *      -Wconversion option set since the sizeof(int) == 8 case will give     *
 *      code that is for a long int on most systems, and hence there will be  *
 *      an implicit conversion. To avoid this, we use the macros defined in   *
 *      limits.h so that the compiler only sees code with legal values for    *
 *      the various data types.                                               *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_bytes.h:                                                 *
 *          Header where rssringoccs_Endian enum data type is defined, and    *
 *          where the function's prototype is given.                          *
 *  2.) limits.h:                                                             *
 *          Standard C library header file containing the CHAR_BIT macro      *
 *          which tells us how many bits are in one char.                     *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in config_librssringoccs.sh uses gcc and has the  *
 *  -pedantic and -std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 14, 2021                                              *
 ******************************************************************************/

/*  Standard library file containing the CHAR_BIT macro and more.             */
#include <limits.h>

/*  Where the function's prototype and rssringoccs_Endian are defined.        */
#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

/******************************************************************************
 *                             The Basic Idea                                 *
 ******************************************************************************
 *  We'll use union to have an unsigned int and a char array share the same   *
 *  memory. We'll then set the integer part to the hexidecimal number         *
 *  0x01020304 (for platforms with sizeof(int) = 4). The char array will      *
 *  see this as follows for little-endian systems:                            *
 *                                                                            *
 *          -----------------------------                                     *
 *          |  04  |  03  |  02  |  01  |                                     *
 *          -----------------------------                                     *
 *                                                                            *
 *  Whereas on big-endian the char array will see:                            *
 *                                                                            *
 *          -----------------------------                                     *
 *          |  01  |  02  |  03  |  04  |                                     *
 *          -----------------------------                                     *
 *                                                                            *
 *  By checking the zeroth element of the array, we can determine endianness. *
 ******************************************************************************/

/*  Check how many bits in a char. This is usually 8. The C89/C90 standard,   *
 *  as well as the C99 and C11 standards, require a char to be at least 8     *
 *  bits, but the data type can be more. rss_ringoccs only works for systems  *
 *  supporting 8 and 16 bit chars. Anything else will cause a compiler error. */
#if CHAR_BIT == 8

/*  Function for determining the endianness of the system.                    */
rssringoccs_Endian rssringoccs_Determine_Endianness(void)
{

    /*  For 16-bit integer, the maximum value is 32767. Check this. On most   *
     *  platforms (all that I've tested) 16-bit is usually a short, whereas   *
     *  an int is 32-bit. The C89/C90 standard only requires int to be at     *
     *  least 16-bit, so this scenario is "technically" possible.             */
#if INT_MAX == 32767

    /*  Use the union C-keyword to create a data-type that has an unsigned    *
     *  int and a char array consisting of two elements, which share the same *
     *  address in memory. By setting the integer portion to the hexidecimal  *
     *  0x0102 we can use the fact that the char array c is occupying the     *
     *  same address as the unsigned int i and interpret the value as a char. *
     *  With this we can see if the zeroth value of the array is 01 or if it  *
     *  is 02. This will tell us if we have little-endian or big-endian.      */
    union {
        unsigned int i;
        char c[2];
    } e = { 0x0102 };

    /*  Check the value of the zeroth element of the char array. If it is 1,  *
     *  then we know we have big-endian. If it is 2, we have little-endian.   */
    if (e.c[0] == 1)
        return rssringoccs_BigEndian;
    else
        return rssringoccs_LittleEndian;

    /*  A 32-bit signed int has maximum value 2147483647. Check this. This is *
     *  the usual size encountered for int.                                   */
#elif INT_MAX == 2147483647

    /*  Same as before, use union to define an unsigned int, and trick the    *
     *  computer into thinking the address has the values for a char array.   *
     *  Read of the zeroth entry to determine endianness.                     */
    union {
        unsigned int i;
        char c[4];
    } e = { 0x01020304 };

    /*  If the zeroth entry is 2 or 3, we have mixed-endian. If 1, we have    *
     *  big-endian, and lastly 4 means little-endian.                         */
    switch (e.c[0])
    {
        case 0x01:
            return rssringoccs_BigEndian;
        case 0x04:
            return rssringoccs_LittleEndian;
        default:
            return rssringoccs_MixedEndian;
    }
    
    /*  The last size we'll check is 64-bit integer types. This is almost     *
     *  certainly not the case since long and long long (a C99 data-type)     *
     *  are usually 64-bit, but it never hurts to have the code available.    */
#elif INT_MAX == 9223372036854775807

    /*  Define a union of a char array and an unsigned int and use this       *
     *  to determine endianness.                                              */
    union {
        unsigned int i;
        char c[8];
    } e = { 0x0102030405060708 };

    /*  Big-endian will correspond to 1, little-endian to 8, and mixed        *
     *  will be a value between 2 and 7, inclusive.                           */
    switch (e.c[0])
    {
        case 0x01:
            return rssringoccs_BigEndian;
        case 0x08:
            return rssringoccs_LittleEndian;
        default:
            return rssringoccs_MixedEndian;
    }
    
    /*  And here we give up and return unknown. This breaks portability, but  *
     *  only for very bizarre systems. If your system gets here, feel free to *
     *  ask us to add in the details or edit the code yourself. There may     *
     *  exist the strange 24-bit, 36-bit, and 48-bit ints, who knows.         */
#else
    return rssringoccs_UnknownEndian;
#endif
    /*  End of #if INT_MAX == 32767.                                          */
}
/*  End of rssringoccs_Determine_Endianness definition with CHAR_BIT == 8.    */

/*  We'll allow for 16-bit chars as well, since they seem to exist. Any other *
 *  size we'll simply raise error. If your platform is in this bizarre case   *
 *  feel free to edit the code. It's all licensed under GPL 3.                */
#elif CHAR_BIT == 16
/*  Else statement of #if CHAR_BIT == 8.                                      */

/*  Function for determining the endianness of the system.                    */
rssringoccs_Endian rssringoccs_Determine_Endianness(void)
{
    /*  If int and char are the same size, our trick won't work. Return       *
     *  unknown.                                                              */
#if INT_MAX == 32767
    return rssringoccs_UnknownEndian;

    /*  If we have a 32-bit int with 16-bit char, we apply the same trick as  *
     *  before using a char array with two elements.                          */
#elif INT_MAX == 2147483647
    union {
        unsigned int i;
        char c[2];
    } e = { 0x00010002 };

    if (e.c[0] == 1)
        return rssringoccs_BigEndian;
    else
        return rssringoccs_LittleEndian;

    /*  Same idea but for 64-bit int and 16-bit char.                         */
#elif INT_MAX == 9223372036854775807
    union {
        unsigned int i;
        char c[4];
    } e = { 0x0001000200030004 };

    switch (e.c[0])
    {
        case 0x0001:
            return rssringoccs_BigEndian;
        case 0x0004:
            return rssringoccs_LittleEndian;
        default:
            return rssringoccs_MixedEndian;
    }

    /*  Finally, 128-bit int with 16-bit char, which I'm pretty sure          *
     *  doesn't exist, but why not, here's the code.                          */
#elif INT_MAX == 18446744073709551615
    union {
        unsigned int i;
        char c[4];
    } e = { 0x00010002000300040005000600070008 };

    switch (e.c[0])
    {
        case 0x0001:
            return rssringoccs_BigEndian;
        case 0x0008:
            return rssringoccs_LittleEndian;
        default:
            return rssringoccs_MixedEndian;
    }
    
    /*  Anything else, return unknown.                                        */
#else
    return rssringoccs_UnknownEndian;
#endif
    /*  End of #if INT_MAX == 32767.                                          */
}
/*  End of rssringoccs_Determine_Endianness with CHAR_BIT == 16.              */

/*  For CHAR_BIT not equal to 8 of 16 abort the compiling. rss_ringoccs does  *
 *  not support such platforms. Only the RSR reader needs this code, so the   *
 *  rest of rss_ringoccs can be used. Just remove this from the build.        */
#else
/*  Else statement for #if CHAR_BIT == 8.                                     */

#error "CHAR_BIT is neither 8 nor 16. Aborting compilation."

#endif
/*  End of #if CHAR_BIT == 8.                                                 */
