/*  Main header for the math functions. Contains <math.h> as well.            */
#include "special_functions.h"

float Normeq_Float(float *w_func, long n_elements){
    float tot_sq = 0.0;
    float sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += w_func[i];
        sq_tot += w_func[i]*w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

double Normeq_Double(double *w_func, long n_elements){
    double tot_sq = 0.0;
    double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += w_func[i];
        sq_tot += w_func[i]*w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

long double Normeq_Long_Double(long double *w_func, long n_elements){
    long double tot_sq = 0.0;
    long double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += w_func[i];
        sq_tot += w_func[i]*w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

double Normeq_Short(short *w_func, long n_elements){
    double tot_sq = 0.0;
    double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += (double)w_func[i];
        sq_tot += (double)w_func[i] * (double)w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

double Normeq_Int(int *w_func, long n_elements){
    double tot_sq = 0.0;
    double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += (double)w_func[i];
        sq_tot += (double)w_func[i] * (double)w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

double Normeq_Long(long *w_func, long n_elements){
    double tot_sq = 0.0;
    double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += (double)w_func[i];
        sq_tot += (double)w_func[i] * (double)w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}

double Normeq_Long_Long(long long *w_func, long n_elements){
    double tot_sq = 0.0;
    double sq_tot = 0.0;
    long i;

    for (i = 0; i < n_elements; i++) {

        /* Compute the sum and the sum of the squares.                        */
        tot_sq += (double)w_func[i];
        sq_tot += (double)w_func[i] * (double)w_func[i];
    }

    /*  Square the sum.                                                       */
    tot_sq *= tot_sq;
    return n_elements * sq_tot / tot_sq;
}
