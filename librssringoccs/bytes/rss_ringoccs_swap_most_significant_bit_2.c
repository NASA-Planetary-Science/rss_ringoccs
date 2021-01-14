
#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

void rssringoccs_Swap_Most_Significant_Bit_2(char *ptr)
{
    rssringoccs_Swap_Bytes(&(ptr[0]), &(ptr[1]));
}
