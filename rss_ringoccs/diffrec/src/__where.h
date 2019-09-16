#ifndef RSS_RINGOCCS_WHERE_H
#define RSS_RINGOCCS_WHERE_H
#include <stdlib.h>

extern long **Where_Greater_Char(char *data, long dim, double threshold);

extern long **
Where_Greater_UChar(unsigned char *data, long dim, double threshold);

extern long **Where_Greater_Short(short *data, long dim, double threshold);

extern long **
Where_Greater_UShort(unsigned short *data, long dim, double threshold);

extern long **Where_Greater_Int(int *data, long dim, double threshold);

extern long **
Where_Greater_UInt(unsigned int *data, long dim, double threshold);

extern long **Where_Greater_Long(long *data, long dim, double threshold);

extern long **
Where_Greater_ULong(unsigned long *data, long dim, double threshold);

extern long **
Where_Greater_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Greater_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **Where_Greater_Float(float *data, long dim, double threshold);

extern long **Where_Greater_Double(double *data, long dim, double threshold);

extern long **
Where_Greater_Long_Double(long double *data, long dim, long double threshold);

extern long **Where_Lesser_Char(char *data, long dim, double threshold);

extern long **
Where_Lesser_UChar(unsigned char *data, long dim, double threshold);

extern long **Where_Lesser_Short(short *data, long dim, double threshold);

extern long **
Where_Lesser_UShort(unsigned short *data, long dim, double threshold);

extern long **Where_Lesser_Int(int *data, long dim, double threshold);

extern long **
Where_Lesser_UInt(unsigned int *data, long dim, double threshold);

extern long **Where_Lesser_Long(long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong(unsigned long *data, long dim, double threshold);

extern long **Where_Lesser_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **Where_Lesser_Float(float *data, long dim, double threshold);

extern long **Where_Lesser_Double(double *data, long dim, double threshold);

extern long **
Where_Lesser_Long_Double(long double *data, long dim, long double threshold);

#endif