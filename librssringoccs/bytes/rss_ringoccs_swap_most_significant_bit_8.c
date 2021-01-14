
#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

void rssringoccs_Swap_Most_Significant_Bit_8(char *ptr)
{
    rssringoccs_Swap_Bytes(&(ptr[0]), &(ptr[7]));
    rssringoccs_Swap_Bytes(&(ptr[1]), &(ptr[6]));
    rssringoccs_Swap_Bytes(&(ptr[2]), &(ptr[5]));
    rssringoccs_Swap_Bytes(&(ptr[3]), &(ptr[4]));
}
