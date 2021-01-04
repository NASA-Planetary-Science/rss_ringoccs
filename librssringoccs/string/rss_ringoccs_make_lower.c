
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>

void rssringoccs_Make_Lower(char *str)
{
    unsigned long n = 0;

    if (str == NULL)
        return;

    while(str[n])
    {
        str[n] = tolower(str[n]);
        ++n;
    }
}