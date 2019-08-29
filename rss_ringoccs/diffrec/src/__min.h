#ifndef RSS_RINGOCCS_MIN_H
#define RSS_RINGOCCS_MIN_H

float Min_Float(float *arr, long n_elements){
    float min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

double Min_Double(double *arr, long n_elements){
    double min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

long double Min_Long_Double(long double *arr, long n_elements){
    long double min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

short Min_Short(short *arr, long n_elements){
    short min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

int Min_Int(int *arr, long n_elements){
    int min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

long Min_Long(long *arr, long n_elements){
    long min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

long Min_Long_Long(long long *arr, long n_elements){
    long long min = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}

#endif