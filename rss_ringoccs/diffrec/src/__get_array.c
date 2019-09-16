#include "__get_array.h"

void Get_Char_to_Double_Array(char *x, double *y, long dim, double (*f)(char))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_UChar_to_Double_Array(unsigned char *x, double *y,
                               long dim, double (*f)(unsigned char))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Short_to_Double_Array(short *x, double *y,
                               long dim, double (*f)(short))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_UShort_to_Double_Array(unsigned short *x, double *y,
                                long dim, double (*f)(unsigned short))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Int_to_Double_Array(int *x, double *y, long dim, double (*f)(int))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_UInt_to_Double_Array(unsigned int *x, double *y,
                              long dim, double (*f)(unsigned int))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Long_to_Double_Array(long *x, double *y, long dim, double (*f)(long))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_ULong_to_Double_Array(unsigned long *x, double *y,
                               long dim, double (*f)(unsigned long))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Long_Long_to_Double_Array(long long *x, double *y,
                                   long dim, double (*f)(long long))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_ULong_Long_to_Double_Array(unsigned long long *x, double *y,
                                    long dim, double (*f)(unsigned long long))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Float_Array(float *x, float *y, long dim, float (*f)(float))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Double_Array(double *x, double *y, long dim, double (*f)(double))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Long_Double_Array(long double *x, long double *y, long dim,
                           long double (*f)(long double))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}