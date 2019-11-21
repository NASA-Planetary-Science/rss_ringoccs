#ifndef RSS_RINGOCCS_GET_ARRAY_H
#define RSS_RINGOCCS_GET_ARRAY_H

#include <complex.h>

extern void Get_Char_to_Double_Array(char *x, double *y,
                                     long dim, double (*f)(char));

extern void Get_UChar_to_Double_Array(unsigned char *x, double *y,
                                      long dim, double (*f)(unsigned char));

extern void Get_Short_to_Double_Array(short *x, double *y,
                                      long dim, double (*f)(short));

extern void Get_UShort_to_Double_Array(unsigned short *x, double *y,
                                       long dim, double (*f)(unsigned short));

extern void Get_Int_to_Double_Array(int *x, double *y, long dim,
                                    double (*f)(int));

extern void Get_UInt_to_Double_Array(unsigned int *x, double *y,
                                     long dim, double (*f)(unsigned int));

extern void Get_Long_to_Double_Array(long *x, double *y,
                                     long dim, double (*f)(long));

extern void Get_ULong_to_Double_Array(unsigned long *x, double *y,
                                      long dim, double (*f)(unsigned long));

extern void Get_Long_Long_to_Double_Array(long long *x, double *y,
                                          long dim, double (*f)(long long));

void Get_ULong_Long_to_Double_Array(unsigned long long *x, double *y,
                                    long dim, double (*f)(unsigned long long));

extern void Get_Float_Array(float *x, float *y,
                            long n_elemnts, float (*f)(float));

extern void Get_Double_Array(double *x, double *y, long n_elements,
                             double (*f)(double));

extern void Get_Long_Double_Array(long double *x, long double *y,
                                  long n_elements,
                                  long double (*f)(long double));

#endif