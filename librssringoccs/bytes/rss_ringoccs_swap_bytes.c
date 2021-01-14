
#include <rss_ringoccs/include/rss_ringoccs_bytes.h>

void rssringoccs_Swap_Bytes(char *ptr1, char *ptr2)
{
    char tmp;

    tmp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = tmp;
    return;
}
