#ifndef RSS_RINGOCCS_WINDOW_NORMALIZATION_H
#define RSS_RINGOCCS_WINDOW_NORMALIZATION_H

#include "__math_constants.h"
#include <math.h>

float Window_Normalization_Float(float *ker, long dim, float dx, float f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    float T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += ker[i];
    }
    T1 = fabsf(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

double Window_Normalization_Double(double *ker, long dim,
                                   double dx, double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += ker[i];
    }
    T1 = fabs(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

long double Window_Normalization_Long_Double(long double *ker, long dim,
                                             long double dx,
                                             long double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    long double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += ker[i];
    }
    T1 = fabsl(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

double Window_Normalization_Short(short *ker, long dim,
                                  double dx, double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += (double)ker[i];
    }
    T1 = fabs(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

double Window_Normalization_Int(int *ker, long dim, double dx, double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += (double)ker[i];
    }
    T1 = fabs(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

double Window_Normalization_Long(long *ker, long dim,
                                 double dx, double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += (double)ker[i];
    }
    T1 = fabs(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

double Window_Normalization_Long_Long(long long *ker, long dim,
                                      double dx, double f_scale){

    /*  Declare variable for indexing.                                        */
    long i;

    /*  Compute the Free-Space integral.                                      */
    double T1 = 0.0;

    for (i=0; i<dim; ++i){
        T1 += (double)ker[i];
    }
    T1 = fabs(T1 * dx);

    /* Retur the normalization factor.                                        */
    return SQRT_2 * f_scale / T1;
}

#endif