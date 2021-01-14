#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

#include <limits.h>

#if CHAR_BIT == 8

rssringoccs_Endian rssringoccs_Determine_Endianness(void)
{
    if (sizeof(unsigned int) == 2)
    {
        union {
            unsigned int i;
            char c[2];
        } e = { 0x0102 };

        if (e.c[0] == 1)
            return rssringoccs_BigEndian;
        else
            return rssringoccs_LittleEndian;
    }
    else if (sizeof(unsigned int) == 4)
    {
        union {
            unsigned int i;
            char c[4];
        } e = { 0x01020304 };
        switch (e.c[0])
        {
            case 0x01:
                return rssringoccs_BigEndian;
            case 0x04:
                return rssringoccs_LittleEndian;
            default:
                return rssringoccs_MixedEndian;
        }
    }
    else if (sizeof(unsigned int) == 8)
    {
        union {
            unsigned int i;
            char c[8];
        } e = { 0x0102030405060708 };

        switch (e.c[0])
        {
            case 0x01:
                return rssringoccs_BigEndian;
            case 0x08:
                return rssringoccs_LittleEndian;
            default:
                return rssringoccs_MixedEndian;
        }
    }
    else
        return rssringoccs_UnknownEndian;
}

#elif CHAR_BIT == 16

{
    if (sizeof(unsigned int) == 2)
    {
        union {
            unsigned int i;
            char c[2];
        } e = { 0x00010002 };
        if (e.c[0] == 1)
            return rssringoccs_BigEndian;
        else
            return rssringoccs_LittleEndian;
    }
    else if (sizeof(unsigned int) == 4)
    {
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
    }
    else if (sizeof(unsigned int) == 8)
    {
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
    }
    else
        return rssringoccs_UnknownEndian;
}

#else

#error "CHAR_BIT is neither 8 nor 16. Aborting compilation."

#endif
