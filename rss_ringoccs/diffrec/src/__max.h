#ifndef RSS_RINGOCCS_MAX_H
#define RSS_RINGOCCS_MAX_H

float Max_Float(float *arr, long n_elements){
    float max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

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

double Max_Long_Double(long double *arr, long n_elements){
    long double max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

short Max_Short(short *arr, long n_elements){
    short max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

int Max_Int(int *arr, long n_elements){
    int max = arr[0];
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

long long Max_Long_Long(long long *arr, long n_elements){
    long long max = arr[0];
    long i;
    for (i = 1; i < n_elements; ++i){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

#endif