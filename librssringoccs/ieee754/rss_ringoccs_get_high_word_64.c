#include <rss_ringoccs/include/rss_ringoccs_ieee754.h>

unsigned long rssringoccs_Get_High_Word64(rssringoccs_IEE754_Word64 x)
{
	unsigned long out = x.integer;
	out = out >> 52;
	out = out & 2047UL;

	return out;
}
