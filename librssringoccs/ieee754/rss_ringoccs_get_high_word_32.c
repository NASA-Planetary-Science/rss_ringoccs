#include <rss_ringoccs/include/rss_ringoccs_ieee754.h>

unsigned int rssringoccs_Get_High_Word32(rssringoccs_IEE754_Word32 x)
{
	unsigned int out = x.integer;
	out = out >> 23;
	out = out & 255U;

	return out;
}
