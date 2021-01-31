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
 *                           rss_ringoccs_ieee754                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains tools for manipulating floating point variables for          *
 *      platforms using the IEEE754 format (most do in 2021). All of this is  *
 *      experimental and mostly for the sake of learning, and is not directly *
 *      required for rss_ringoccs. Uses of this header file include accurate  *
 *      and fast logarithmic functions, square roots, and more.               *
 *  NOTES:                                                                    *
 *      While the code is written in ANSI C, this is NOT portable since it    *
 *      assumes various things. This part of the code makes the following     *
 *      assumptions (however, as stated, this code is not used for the        *
 *      primary diffraction correction routines):                             *
 *          1.) Your platform uses IEEE754 format for floating point          *
 *              arithmetic. Most modern computers do.                         *
 *          2.) You have 32-bit float and 64-bit double. This is NOT          *
 *              required by the C89/C90 standard, only minimum sizes are      *
 *              specified. 32-bit/64-bit single and double precision is the   *
 *              most common, but this can break portability.                  *
 *          2.) And unsigned long int has 64-bits, and an unsigned int has    *
 *              32-bits. This will most likely NOT be true on 32-bit systems, *
 *              especially 32-bit Microsoft Windows machines where unsigned   *
 *              long int is 32-bit. This assumption is true for 64-bit        *
 *              computers, including all of those rss_ringoccs was tested on. *
 *      Endianness shouldn't matter, however the code has only been tested on *
 *      Little Endian systems.                                                *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 22, 2021                                              *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_IEEE754_H__
#define __RSS_RINGOCCS_IEEE754_H__

/*  Data type for a 32-bit floating point number. This is assumed to          *
 *  correspond to the float data type. Note that char is assumed to be 8 bits.*
 *  Also, an unsigned int is assumed to be 32-bit.                            */
typedef union word32 {
	float real;
	unsigned int integer;
	unsigned char bits[4];
} rssringoccs_IEE754_Word32;

/*  Data type for a 64-bit floating point number. This is assumed to          *
 *  correspond to the double data type. Note that char is assumed to be 8     *
 *  bits. Also, an unsigned long is assumed to be 32-bit.                     */
typedef union word64 {
	double real;
	unsigned long integer;
	unsigned char bits[8];
} rssringoccs_IEE754_Word64;

unsigned int rssringoccs_Get_High_Word32(rssringoccs_IEE754_Word32 x);
unsigned int rssringoccs_Get_Low_Word32(rssringoccs_IEE754_Word32 x);

unsigned long rssringoccs_Get_High_Word64(rssringoccs_IEE754_Word64 x);
unsigned long rssringoccs_Get_Low_Word64(rssringoccs_IEE754_Word64 x);

#endif
/*  End of #ifndef __RSS_RINGOCCS_IEEE754_H__.                                */
