#ifndef RSS_RINGOCCS_MAX_H
#define RSS_RINGOCCS_MAX_H

double Max_Double(double *arr, long n_elements){
    double max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

long Max_Long(long *arr, long n_elements){
    long max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

#endif