#ifndef RSS_RINGOCCS_GET_ARRAY_H
#define RSS_RINGOCCS_GET_ARRAY_H

extern void Get_Float_Array(float *x, float *y,
                            long n_elemnts, float (*f)(float));

extern void Get_Double_Array(double *x, double *y, long n_elements,
                             double (*f)(double));

extern void Get_Long_Double_Array(long double *x, long double *y,
                                  long n_elements,
                                  long double (*f)(long double));

#endif