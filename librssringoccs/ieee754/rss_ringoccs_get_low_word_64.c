#include <rss_ringoccs/include/rss_ringoccs_ieee754.h>

unsigned long rssringoccs_Get_Low_Word64(rssringoccs_IEE754_Word64 x)
{
	unsigned long out = x.integer;
	out = out & 4503599627370495UL;
	return out;
}
