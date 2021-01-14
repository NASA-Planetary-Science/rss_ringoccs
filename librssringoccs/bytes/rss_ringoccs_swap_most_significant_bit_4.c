
#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

void rssringoccs_Swap_Most_Significant_Bit_4(char *ptr)
{
    rssringoccs_Swap_Bytes(&(ptr[0]), &(ptr[3]));
    rssringoccs_Swap_Bytes(&(ptr[1]), &(ptr[2]));
}
