
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Standard algorithm of time-complexity O(n) to take an array and reverse   *
 *  the order. This is equivalent to taking a numpy array arr and writing     *
 *  arr = arr[::-1] (in Python 3, at least). Since all of the main routines   *
 *  use double pointers, only a double version is provided.                   */
void ReverseDoubleArray(double *arr, long arrsize)
{
    double val;
    long i;
    for(i=0; i<arrsize/2; i++)
    {
        val = arr[i];
        arr[i] = arr[arrsize-i-1];
        arr[arrsize-i-1] = val;
    }
}
