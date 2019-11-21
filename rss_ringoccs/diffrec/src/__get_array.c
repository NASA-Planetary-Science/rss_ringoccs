#include "__get_array.h"

/*  Make sure the name __get_where_pointer is available.                      */
#ifdef __get_one_real_from_one_real
#undef __get_one_real_from_one_real
#endif

/*  To avoid repeating the same code over and over again, define this macro   *
 *  to be used for all of the where_lesser functions. Since the only thing    *
 *  that changes between the various functions is the type of the input       *
 *  pointer, the code is exactly the same.                                    */

#define __get_one_real_from_one_real(x, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i]);\
    }\
})

void Get_Char_to_Double_Array(char *x, double *y, long dim, double (*f)(char))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_UChar_to_Double_Array(unsigned char *x, double *y,
                               long dim, double (*f)(unsigned char))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Short_to_Double_Array(short *x, double *y,
                               long dim, double (*f)(short))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_UShort_to_Double_Array(unsigned short *x, double *y,
                                long dim, double (*f)(unsigned short))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Int_to_Double_Array(int *x, double *y, long dim, double (*f)(int))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_UInt_to_Double_Array(unsigned int *x, double *y,
                              long dim, double (*f)(unsigned int))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Long_to_Double_Array(long *x, double *y, long dim, double (*f)(long))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_ULong_to_Double_Array(unsigned long *x, double *y,
                               long dim, double (*f)(unsigned long))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Long_Long_to_Double_Array(long long *x, double *y,
                                   long dim, double (*f)(long long))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_ULong_Long_to_Double_Array(unsigned long long *x, double *y,
                                    long dim, double (*f)(unsigned long long))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Float_Array(float *x, float *y, long dim, float (*f)(float))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Double_Array(double *x, double *y, long dim, double (*f)(double))
{
    __get_one_real_from_one_real(x, y, dim, f);
}

void Get_Long_Double_Array(long double *x, long double *y, long dim,
                           long double (*f)(long double))
{
    __get_one_real_from_one_real(x, y, dim, f);
}
