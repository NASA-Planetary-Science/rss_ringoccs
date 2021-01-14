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
 *                           rss_ringoccs_bytes                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide tools for swapping bytes of char pointers and for determining *
 *      the endianness of your platform. This is used with the RSR reader.    *
 *      The RSR files are big endian whereas most processor architectures use *
 *      little-endian. Indeed, all computers rss_ringoccs were little-endian. *
 *      The tools here allow one to read the RSR files regardless.            *
 *  NOTES:                                                                    *
 *      It is assumed a char is either 8 bits or 16 bits. This is probably    *
 *      the case for your computer. The only devices I know of where this may *
 *      not be true are hand-held calculators. So don't try to use            *
 *      rss_ringoccs on those, I suppose. If CHAR_BITS (defined in limits.h)  *
 *      is neither of these, rss_ringoccs will fail to build.                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 14, 2020                                              *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_BYTES_H__
#define __RSS_RINGOCCS_BYTES_H__

/*  Data type for determining the endianness of your platform. If mixed or    *
 *  unknown endianness are returned by rssringoccs_Determine_Endianness,      *
 *  rss_ringoccs will be unable to extract the data from the RSR files.       */
typedef enum {
    rssringoccs_UnknownEndian,
    rssringoccs_LittleEndian,
    rssringoccs_MixedEndian,
    rssringoccs_BigEndian
} rssringoccs_Endian;

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Determine_Endianness                                      *
 *  Purpose:                                                                  *
 *      Determines the endianness of your platform.                           *
 *  Arguments:                                                                *
 *      None (void).                                                          *
 *  Output:                                                                   *
 *      rssringoccs_Endian endianness:                                        *
 *          A rssringoccs_Endian data type whose value corresponds to the     *
 *          endianness of your system.                                        *
 *  NOTES:                                                                    *
 *      1.) This function assumes CHAR_BITS (defined in limits.h) is either   *
 *          8 or 16. If this is not true, rss_ringoccs will fail to build.    *
 *                                                                            *
 *      2.) If mixed endian is determined, rss_ringoccs will be unable to     *
 *          to extract the data from the RSR files. To the best of my         *
 *          knowledge, mixed endian systems are rare/not in use.              *
 ******************************************************************************/
extern rssringoccs_Endian rssringoccs_Determine_Endianness(void);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Swap_Bytes                                                *
 *  Purpose:                                                                  *
 *      Swaps the values of two char pointers.                                *
 *  Arguments:                                                                *
 *      ptr1 (char *):                                                        *
 *          The first pointer to a char.                                      *
 *      ptr2 (char *):                                                        *
 *          The second pointer to a char.                                     *
 *  Output:                                                                   *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Swap_Bytes(char *ptr1, char *ptr2);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Swap_Most_Significant_Bit_2                               *
 *  Purpose:                                                                  *
 *      Changes the endianness of a data type that is two chars long. This    *
 *      is usually a "short". Similarly functions are provided for 4 and 8.   *
 *  Arguments:                                                                *
 *      ptr (char *):                                                         *
 *          A pointer to a char array.                                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Swap_Most_Significant_Bit_2(char *ptr);
extern void rssringoccs_Swap_Most_Significant_Bit_4(char *ptr);
extern void rssringoccs_Swap_Most_Significant_Bit_8(char *ptr);

#endif
/*  End of include guard.                                                     */
