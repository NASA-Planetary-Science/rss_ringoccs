#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    long double max, temp, start;
    long double x, dx;
    long double *y0, *y1;
    int n, N, ind;
    clock_t t1, t2;

    N = 1e8;
    start = 100.0;
    dx = 2.0L*start / N;
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    x = -start;
    t1 = clock();
    for (n=0; n<N; ++n)
    {
        y0[n] = sinl(x);
        x += dx;
    }
    t2 = clock();
    printf("C99: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    x = -start;
    t1 = clock();
    for (n=0; n<N; ++n)
    {
        y1[n] = rssringoccs_LDouble_Sin(x);
        x += dx;
    }
    t2 = clock();
    printf("rss_ringoccs: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    max = 0.0;
    ind = 0;
    for (n=0; n<N; ++n)
    {
        temp = fabsl(y0[n] - y1[n]);
        if (max < temp)
        {
            max = temp;
            ind = n;
        }
    }

    printf("Max Error: %.24Lf\n", max);
    printf("%.16Lf\n%.24Lf\n%.24Lf\n", -6.0L+ind*dx, y0[ind], y1[ind]);
    free(y0);
    free(y1);
    return 0;
}
