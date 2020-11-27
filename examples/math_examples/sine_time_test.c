#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    double max, temp;
    double x, dx;
    double *y0, *y1;
    int n, N;
    clock_t t1, t2;

    N = 1e7;
    dx = TWO_PI / N;
    x = 0.0;
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    t1 = clock();
    for (n=0; n<N; ++n)
    {
        y0[n] = sin(x);
        x += dx;
    }
    t2 = clock();

    printf("C99: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    x = 0.0;
    t1 = clock();
    for (n=0; n<N; ++n)
    {
        y1[n] = rssringoccs_Double_Sin(x);
        x += dx;
    }
    t2 = clock();

    printf("rss_ringoccs: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    max = 0.0;
    for (n=0; n<N; ++n)
    {
        temp = fabs(y0[n] - y1[n]);
        if (max < temp)
            max = temp;
    }

    printf("Max Error: %.16f\n", max);

    free(y0);
    free(y1);
    return 0;
}
