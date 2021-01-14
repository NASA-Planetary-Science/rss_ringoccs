

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_BYTES_H__
#define __RSS_RINGOCCS_BYTES_H__

/*  We prepend rssringoccs onto False, True, and Bool to avoid name conflicts.*/
typedef enum {
    rssringoccs_UnknownEndian,
    rssringoccs_LittleEndian,
    rssringoccs_MixedEndian,
    rssringoccs_BigEndian
} rssringoccs_Endian;

extern rssringoccs_Endian rssringoccs_Determine_Endianness(void);

extern void rssringoccs_Swap_Bytes(char *ptr1, char *ptr2);

extern void rssringoccs_Swap_Most_Significant_Bit_8(char *ptr);
extern void rssringoccs_Swap_Most_Significant_Bit_4(char *ptr);
extern void rssringoccs_Swap_Most_Significant_Bit_2(char *ptr);

#endif
/*  End of include guard.                                                     */
