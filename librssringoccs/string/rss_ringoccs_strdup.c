

#include <stdlib.h>
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>

char *rssringoccs_strdup(const char *str)
{
    /*  Create space for the output. The +1 is for the NULL terminator.       */
    char *dst = malloc(strlen(str) + 1);

    /*  Check if malloc failed.                                               */
    if (dst == NULL)
        return NULL;

    /*  Copy the input string into dst.                                       */
    strcpy(dst, str);
    return dst;
}
