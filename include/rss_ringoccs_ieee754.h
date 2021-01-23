





typedef union word64 {
	double real;
	unsigned long integer;
	unsigned char bits[8];
} rssringoccs_IEE754_Word64;

unsigned long rssringoccs_Get_High_Word64(rssringoccs_IEE754_Word64 x);

unsigned long rssringoccs_Get_Low_Word64(rssringoccs_IEE754_Word64 x);