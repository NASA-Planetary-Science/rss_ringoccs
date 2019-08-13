#ifndef RSS_RINGOCCS_MATH_FUNCTIONS_H
#define RSS_RINGOCCS_MATH_FUNCTIONS_H

double max_double(double *arr, long n_elements){
    double max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

long max_long(long *arr, long n_elements){
    long max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

double min_double(double *arr, long n_elements){
    double min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

long min_long(long *arr, long n_elements){
    long min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

#endif