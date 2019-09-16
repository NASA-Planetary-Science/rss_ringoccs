#include "__window_width.h"

long **WhereGreater(double *data, long dim, double threshold)
{
    long i;
    long arr_size = 0;
    long *dummy = (long *)malloc(sizeof(long) * dim);
    for (i=0; i<dim; ++i) {
        if (data[i] > threshold){
            dummy[arr_size] = i;
            arr_size += 1;
        }
    }
    long *arr = (long *)malloc(sizeof(long) * arr_size);
    for (i=0; i<arr_size; ++i){
        arr[i] = dummy[i];
    }
    free(dummy);
    long **where = (long **)malloc(sizeof(long) * 2);
    where[0] = arr;
    where[1] = &arr_size;
    return where;
}